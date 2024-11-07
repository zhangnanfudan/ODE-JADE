# R Scripts

This directory contains R scripts that define functions and utilities used in the experiments.

## Files

- **data_generate_additive_scaled.R**: Generates data for additive models with scaling.
- **bcd_additive.R**: Implements the BCD algorithm for additive models.
- **twostep_additive.R**: Contains the two-step additive method for model fitting.
- **twostep_additive_init.R**: Initializes parameters for the two-step additive method.

## Usage

These scripts are sourced by the main experiment scripts in the `experiments` directory. They define functions and methods that are called during the execution of experiments.

## Dependencies

Ensure you have the necessary R packages installed, such as `fda` and `Matrix`.
