# SymDock Funnel Plotting Script

**Author:** Allon Goldberg  
**Position:** Research Assistant, Flatiron Institute  
**Date:** May 2025

This repository contains a Python script for generating funnel plots from [SymDock](https://www.rosettacommons.org/docs/latest/application_documentation/docking/SymmetricDocking) simulations. The plots compare global and local interface refinement, with each funnel colored by Fnat to assess docking quality.

## Purpose

The goal of this script is to evaluate docking performance by visualizing energy funnels (score vs. RMSD), separated by interface refinement strategy:

- Global refinement
- Local refinement

Each point in the plot is colored according to its Fnat value (fraction of native contacts recovered), providing a qualitative view of docking accuracy.

## Usage

This script is executed from the command line and takes a single string argument. The argument is not used internally and serves only to satisfy a required input format.

## Required Directory Structure

The script assumes a specific folder structure for input score files. Each numbered directory (e.g., `3/`, `5/`) represents a separate SymDock run and must contain both global and local refinement results:

```
RUN_DIRECTORY/
├── 3/
│   ├── score.sc
│   └── LOCAL/
│   └── score.sc
└── 5/
   ├── score.sc
   └── LOCAL/
   └── score.sc
```

## Plot Details

- X-axis: RMSD
- Y-axis: Rosetta total score
- Color: Fnat (fraction of native contacts)

The resulting plots help identify favorable docking models and distinguish interface quality between global and local refinement.

## Dependencies

This script requires the following Python packages:

- `matplotlib`
- `pandas`
- `numpy`
- `seaborn` (optional, for enhanced styling)

Install them using pip:
