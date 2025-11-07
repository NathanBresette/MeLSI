# MeLSI Results Reproducibility Scripts

This folder contains R scripts to reproduce each results table from the MeLSI paper.

## Prerequisites

Install required R packages:

```r
install.packages(c("vegan", "ggplot2", "dplyr", "GUniFrac"))

# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("microbiome", "phyloseq"))

# Install MeLSI package
devtools::install_github("NathanBresette/MeLSI")
```

## Scripts

| Script | Table | Description |
|--------|-------|-------------|
| `table1_type1_error.R` | Table 1 | Type I Error Control on Null Data |
| `table2_power_analysis.R` | Table 2 | Method Comparison on Synthetic and Real Datasets |
| `table2_dietswap.R` | (supplement) | DietSwap-only reproduction of Table 2 real-data row |
| `table3_scalability.R` | Table 3 | Scalability Across Sample Size and Dimensionality |
| `table4_parameter_sensitivity.R` | Table 4 | Parameter Sensitivity Analysis |
| `table5_prefiltering.R` | Table 5 | Benefit of Conservative Pre-filtering |

## Usage

Each script is standalone and can be run independently:

```bash
Rscript table1_type1_error.R
Rscript table2_power_analysis.R
# ... etc
```

Or run all scripts at once:

```bash
Rscript run_all_tables.R
```

## Expected Runtime

- Table 1: ~5 minutes
- Table 2: ~15-20 minutes (includes synthetic power analysis + real datasets)
- Table 3: ~20-30 minutes (tests multiple dimensions)
- Table 4: ~10-15 minutes (parameter sweeps)
- Table 5: ~10 minutes

## Output

Each script generates:
1. Console output showing progress and results
2. A CSV file with the results (e.g., `table1_results.csv`)
3. Summary statistics and interpretation

## Notes

- All scripts use `set.seed(42)` for reproducibility
- MeLSI uses 200 permutations (n_perms=200, B=30) for stable and accurate p-values
- Comparison methods (Euclidean, Bray-Curtis, etc.) use 999 permutations (standard practice)
- Real data (Atlas1006 and DietSwap) is loaded from the `microbiome` R package
- DietSwap ships with the package; no external downloads are required

