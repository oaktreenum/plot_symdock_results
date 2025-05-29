# SymDock to funnel plot script, 
# Allon Goldberg, Research Assistant, Flatiron Institute, 5/2025

# Usage: python plot_d4.py <any_string>

# NECESSARY FILE STRUCTURE:
# RUN_DIRECTORY/
# ├── 2/
# │   └── score.sc
# │   └── LOCAL/
# │         └── score.sc
# └── 5/
#     └── score.sc
#     └── LOCAL/
#           └── score.sc


import os
import sys
import glob
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.path as mpath
from io import StringIO
import numpy as np
from scipy.special import logsumexp

current_dir = os.path.basename(os.getcwd())

def SymDockData(score_file):
    with open(score_file, 'r') as scores:
        score_lines = scores.readlines()

    p_near = 0
    df_lines = []
    mode = 0
    for line in score_lines:
        # Find where header starts and change mode to read data
        if 'SCORE' in line:
            mode = 1
            df_lines.append(line)
            continue
        # If end of data, change mode
        if 'SCORE' not in line:
            mode = 0
        
    fake_file = '\n'.join(df_lines)
    data_io = StringIO(fake_file)
    df_raw = pd.read_csv(data_io, sep='\s+', index_col=False,on_bad_lines='skip').fillna(0)
    df = df_raw.iloc[:,1:]
    # Minimum Energies for Axes Scaling
    Energy_min = df['total_score'].min()
    RMS_min = df['rms'].min()
    I_sc_min = df['I_sc'].min()
    # WRONG RMSD COLUMN DO NOT USE
    # Sym_RMSD_min = df['symmepentc_rms'].min()

    return df, Energy_min, RMS_min, I_sc_min

def Pnear(df1,df2,score_column,rmsd_column):

    df1 = df1[[score_column,rmsd_column]]
    df2 = df2[[score_column,rmsd_column]]

    # Combine DataFrames for unified processing
    combined_df = pd.concat([df1, df2], ignore_index=True)
    
    # Constants
    lambda_val = 2.5
    kBT = 0.62   

    # Subtract the minimum score for normalization
    combined_df['adjusted_score'] = combined_df[score_column] - combined_df[score_column].min()

    # Compute the numerator terms: -rmsd^2 / lambda^2 - score / kBT
    combined_df['numerator_exp'] = (
        -combined_df[rmsd_column]**2 / lambda_val**2
        - combined_df['adjusted_score'] / kBT
    )

    # Compute the denominator terms: -score / kBT
    combined_df['denominator_exp'] = -combined_df['adjusted_score'] / kBT

    # Use log-sum-exp
    numerator_lse = logsumexp(combined_df['numerator_exp'])
    denominator_lse = logsumexp(combined_df['denominator_exp'])

    # Compute Pnear
    log_pnear = numerator_lse - denominator_lse
    pnear = np.exp(log_pnear)
    return pnear

def main():
    rmsd_cutoff = 20
    df_dim, dim_Energy_min, dim_RMS_min, dim_I_sc_min = SymDockData("2/score.sc")
    dim_zoom = df_dim[(df_dim['symmetric_rms'] < rmsd_cutoff)]
    df_dimL, dim_Energy_minL, dim_RMS_minL, dim_I_sc_minL = SymDockData("2/LOCAL/score.sc")
    dim_zoomL = df_dimL[(df_dimL['symmetric_rms'] < rmsd_cutoff)]
    df_pent, pent_Energy_min, pent_RMS_min, pent_I_sc_min = SymDockData("5/score.sc")
    df_pentL, pentL_Energy_min, pentL_RMS_min, pentL_I_sc_min = SymDockData("5/LOCAL/score.sc")
    pent_zoom = df_pent[(df_pent['symmetric_rms'] < rmsd_cutoff)]
    pent_zoomL = df_pentL[(df_pentL['symmetric_rms'] < rmsd_cutoff)]

    # calculate Pnear
    pnear_dim = Pnear(df_dim,df_dimL,'I_sc','symmetric_rms')
    pnear_pent = Pnear(df_pent,df_pentL,'I_sc','symmetric_rms')

    fig, axes = plt.subplots(1, 2, figsize=(23, 8))
    fig.suptitle(f'{current_dir} SymDock Plots')

    dim_scatter_plot = axes[0].scatter(dim_zoom['symmetric_rms'],dim_zoom['I_sc'], c=dim_zoom['Fnat'], cmap="crest", alpha=.8)
    dimL_scatter_plot = axes[0].scatter(dim_zoomL['symmetric_rms'],dim_zoomL['I_sc'], marker="+", c="black",alpha=.4, s=60)
    axes[0].set_title('Dimer')
    axes[0].set_xlabel('RMSD')
    xticks = np.arange(0, 20.1, 2.5) 
    axes[0].set_xticks(xticks)
    axes[0].set_ylabel('Interface Score')
    axes[0].annotate(f'PNear: {pnear_dim:.4f}', xy=(0.72, 0.043), xycoords='axes fraction', fontsize=13, color='black')
    cbar0 = plt.colorbar(dim_scatter_plot, ax=axes[0])
    cbar0.set_label('Fnat')

    pent_scatter_plot = axes[1].scatter(pent_zoom['symmetric_rms'], pent_zoom['I_sc'], c=pent_zoom['Fnat'], cmap="flare", alpha=.8)
    pentL_scatter_plot = axes[1].scatter(pent_zoomL['symmetric_rms'], pent_zoomL['I_sc'], marker="+", c="black",alpha=.4, s=60)
    axes[1].set_title('Pentamer')
    xticks1 = np.arange(0, 20.1, 2.5) 
    axes[1].set_xticks(xticks1)
    axes[1].set_xlabel('RMSD')
    axes[1].set_ylabel('Interface Score')
    axes[1].annotate(f'PNear: {pnear_pent:.4f}', xy=(0.72, 0.043), xycoords='axes fraction', fontsize=13, color='black')
    cbar1 = plt.colorbar(pent_scatter_plot, ax=axes[1])
    cbar1.set_label('Fnat')

    with open(f"{current_dir}_summary.txt", 'w') as file:
        file.write(f"Dimer Pnear: {pnear_dim}, Pentamer Pnear: {pnear_pent}\n")
        file.write(f"Dimer Energy Min: {dim_Energy_min}, TAG: {df_dim[df_dim['total_score'] == dim_Energy_min]['description'].iloc[0]}\n")
        file.write(f"Dimer RMS Min: {dim_RMS_min}, TAG: {df_dim[df_dim['rms'] == dim_RMS_min]['description'].iloc[0]}\n")
        file.write(f"Dimer I_sc Min: {dim_I_sc_min}, TAG: {df_dim[df_dim['I_sc'] == dim_I_sc_min]['description'].iloc[0]}\n")
        file.write(f"Pentamer RMS Min: {pent_RMS_min}, TAG: {df_pent[df_pent['rms'] == pent_RMS_min]['description'].iloc[0]}\n")
        file.write(f"Pentamer Energy Min: {pent_Energy_min}, TAG: {df_pent[df_pent['total_score'] == pent_Energy_min]['description'].iloc[0]}\n")
        file.write(f"Pentamer I_sc Min: {pent_I_sc_min}, TAG: {df_pent[df_pent['I_sc'] == pent_I_sc_min]['description'].iloc[0]}\n")


    fig.savefig(f'sdplot_{current_dir}.png', dpi=300)
    fig.savefig(f'sdplot_{current_dir}.svg')

    return



if __name__ == "__main__":
    if sys.argv[1] == 'help':
        print("""Usage: python plot_d4.py <any_string> 
                    NECESSARY FILE STRUCTURE:
                    RUN_DIRECTORY/
                    ├── 2/
                    │   └── score.sc
                    │   └── LOCAL/
                    │         └── score.sc
                    └── 5/
                      └── score.sc
                        └── LOCAL/
                            └── score.sc""")
    else:
        logs = sys.argv[1]
        print(f'\n♡♡♡ Thank you for choosing me as your plotter. I work hard to plot your figures ♡♡♡ \n\nPLOTTING {current_dir}...\n')
        main()
        print(f'The plots for {current_dir} were generated and saved as png(300dpi)/svg files. ☺ ENJOY ☺\n')