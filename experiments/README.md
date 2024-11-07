# Running Experiments

This directory contains scripts for running simulation experiments. The main entry point is `main.R`.

## Usage

Run experiments using the command line:

```bash
Rscript main.R <family> <timepoints> <replications> [options]
```

### Required Arguments

1. **Distribution Family** (first argument)
   - `g` for Gaussian
   - `p` for Poisson
   - `b` for Binomial

2. **Number of Time Points** 
   - Format: `N<number>`
   - Example: `N100` means 101 time points
   - Must include the 'N' prefix

3. **Number of Replications**
   - Format: `R<number>`
   - Example: `R5` means 5 replications
   - Must include the 'R' prefix

### Family-Specific Required Arguments

- **Signal-to-Noise Ratio** (Gaussian family only)
  - Format: `SNR<number>` or `SNRInf`
  - Example: `SNR2` for SNR=2, `SNRInf` for infinite SNR
  - Required when family='g'
  - Not used for Poisson or Binomial families

### Optional Arguments

1. **Number of Tolerance Levels**
   - Format: `NT<number>`
   - Default: 1 if not specified
   - Example: `NT3`

2. **Number of Reweighting Steps**
   - Format: `NW<number>`
   - Default: 4 if not specified
   - Example: `NW5`

### Example Commands

1. **Gaussian Distribution**:
   ```bash
   # Basic usage with SNR=2
   Rscript main.R g N100 R5 SNR2

   # With optional parameters
   Rscript main.R g N100 R5 SNR2 NT3 NW5

   # With infinite SNR
   Rscript main.R g N100 R5 SNRInf
   ```

2. **Poisson Distribution**:
   ```bash
   # Basic usage
   Rscript main.R p N50 R3

   # With optional parameters
   Rscript main.R p N50 R3 NT2 NW4
   ```

3. **Binomial Distribution**:
   ```bash
   # Basic usage
   Rscript main.R b N200 R10

   # With optional parameters
   Rscript main.R b N200 R10 NT3 NW6
   ```

### Output

- Results are saved in the `result` directory
- Output filename format: `<family>_a_N<timepoints>R<replications>NT<tolerance>NW<reweight>[SNR<value>].txt`
- For Gaussian family, SNR value is appended to the filename

### Important Notes

1. The order of arguments doesn't matter except for the family parameter (must be first)
2. SNR parameter is required only for Gaussian family
3. All parameters must follow their exact format (e.g., 'N100' not just '100')
4. The script will create the `result` directory if it doesn't exist
5. By default, existing result files will be overwritten

## Error Handling

The script will stop with an error message if:
- Required arguments are missing
- Arguments are in incorrect format
- SNR is missing for Gaussian family
- SNR is provided for non-Gaussian families
