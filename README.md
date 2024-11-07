# Project Overview

This project is designed to run experiments for ``A Joint estimation approach to sparse additive ordinary differential equations''. 

## Directory Structure

- **R/**: R scripts for the main algorithm and auxillary utilities.
- **experiments/**: Scripts for setting up and running simulation studies.
- **application/**: Real data analysis


## Getting Started

### Simulation

1. Navigate to the `experiments` directory.
2. Check `main.R`, modify the script to set your desired experiment parameters.
3. Run `main.R` using `Rscript` with proper arguments in the command line to execute the experiments.
4. Check the `result` directory for output files.

### Real Data Analysis

- Run `yeast_selected.r` or `stock_old_bin3.r`.