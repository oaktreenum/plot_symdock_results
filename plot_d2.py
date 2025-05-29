# SymDock to funnel plot script, 
# Allon Goldberg, Research Assistant, Flatiron Institute, 5/2025

# Usage: python plot_d2.py <any_string>

# NECESSARY FILE STRUCTURE:
# RUN_DIRECTORY/
# ├── 2/
# │   └── score.sc
# │   └── LOCAL/
# │         └── score.sc
# └── 3/
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
    df_tri, tri_Energy_min, tri_RMS_min, tri_I_sc_min = SymDockData("3/score.sc")
    tri_zoom = df_tri[(df_tri['symmetric_rms'] < rmsd_cutoff)]
    df_triL, tri_Energy_minL, tri_RMS_minL, tri_I_sc_minL = SymDockData("3/LOCAL/score.sc")
    tri_zoomL = df_triL[(df_triL['symmetric_rms'] < rmsd_cutoff)]
    df_dim, dim_Energy_min, dim_RMS_min, dim_I_sc_min = SymDockData("2/score.sc")
    df_dimL, dimL_Energy_min, dimL_RMS_min, dimL_I_sc_min = SymDockData("2/LOCAL/score.sc")
    dim_zoom = df_dim[(df_dim['symmetric_rms'] < rmsd_cutoff)]
    dim_zoomL = df_dimL[(df_dimL['symmetric_rms'] < rmsd_cutoff)]

    # calculate Pnear
    pnear_tri = Pnear(df_tri,df_triL,'I_sc','symmetric_rms')
    pnear_dim = Pnear(df_dim,df_dimL,'I_sc','symmetric_rms')

    fig, axes = plt.subplots(1, 2, figsize=(23, 8))
    fig.suptitle(f'{current_dir} SymDock Plots')

    tri_scatter_plot = axes[1].scatter(tri_zoom['symmetric_rms'],tri_zoom['I_sc'], c=tri_zoom['Fnat'], cmap="crest", alpha=.8)
    triL_scatter_plot = axes[1].scatter(tri_zoomL['symmetric_rms'],tri_zoomL['I_sc'], marker="+", c="black",alpha=.4, s=60)
    axes[1].set_title('Trimer')
    axes[1].set_xlabel('RMSD')
    xticks = np.arange(0, 20.1, 2.5) 
    axes[1].set_xticks(xticks)
    axes[1].set_ylabel('Interface Score')
    axes[1].annotate(f'PNear: {pnear_tri:.4f}', xy=(0.72, 0.043), xycoords='axes fraction', fontsize=13, color='black')
    cbar0 = plt.colorbar(tri_scatter_plot, ax=axes[1])
    cbar0.set_label('Fnat')

    dim_scatter_plot = axes[0].scatter(dim_zoom['symmetric_rms'], dim_zoom['I_sc'], c=dim_zoom['Fnat'], cmap="flare", alpha=.8)
    dimL_scatter_plot = axes[0].scatter(dim_zoomL['symmetric_rms'], dim_zoomL['I_sc'], marker="+", c="black",alpha=.4, s=60)
    axes[0].set_title('Dimer')
    xticks1 = np.arange(0, 20.1, 2.5) 
    axes[0].set_xticks(xticks1)
    axes[0].set_xlabel('RMSD')
    axes[0].set_ylabel('Interface Score')
    axes[0].annotate(f'PNear: {pnear_dim:.4f}', xy=(0.72, 0.043), xycoords='axes fraction', fontsize=13, color='black')
    cbar1 = plt.colorbar(dim_scatter_plot, ax=axes[0])
    cbar1.set_label('Fnat')

    with open(f"{current_dir}_summary.txt", 'w') as file:
        file.write(f"Dimer Pnear: {pnear_dim}, Trimer Pnear: {pnear_tri}\n")
        file.write(f"Dimer RMS Min: {dim_RMS_min}, TAG: {df_dim[df_dim['rms'] == dim_RMS_min]['description'].iloc[0]}\n")
        file.write(f"Dimer Energy Min: {dim_Energy_min}, TAG: {df_dim[df_dim['total_score'] == dim_Energy_min]['description'].iloc[0]}\n")
        file.write(f"Dimer I_sc Min: {dim_I_sc_min}, TAG: {df_dim[df_dim['I_sc'] == dim_I_sc_min]['description'].iloc[0]}\n")
        file.write(f"Trimer Energy Min: {tri_Energy_min}, TAG: {df_tri[df_tri['total_score'] == tri_Energy_min]['description'].iloc[0]}\n")
        file.write(f"Trimer RMS Min: {tri_RMS_min}, TAG: {df_tri[df_tri['rms'] == tri_RMS_min]['description'].iloc[0]}\n")
        file.write(f"Trimer I_sc Min: {tri_I_sc_min}, TAG: {df_tri[df_tri['I_sc'] == tri_I_sc_min]['description'].iloc[0]}\n")
        


    fig.savefig(f'sdplot_{current_dir}.png', dpi=300)
    fig.savefig(f'sdplot_{current_dir}.svg')

    return



if __name__ == "__main__":
    if sys.argv[1] == 'help':
        print("""Usage: python plot_d2.py <any_string> 
                    NECESSARY FILE STRUCTURE:
                    RUN_DIRECTORY/
                    ├── 2/
                    │   └── score.sc
                    │   └── LOCAL/
                    │         └── score.sc
                    └── 3/
                      └── score.sc
                        └── LOCAL/
                            └── score.sc""")
    else:
        logs = sys.argv[1]
        print(f'\n♡♡♡ Thank you for choosing me as your plotter. I work hard to plot your figures ♡♡♡ \n\nPLOTTING {current_dir}...\n')
        main()
        print(f'The plots for {current_dir} were generated and saved as png(300dpi)/svg files. ☺ ENJOY ☺\n')