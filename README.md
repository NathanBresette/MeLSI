# MeLSI: Metric Learning for Statistical Inference in Microbiome Analysis

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R](https://img.shields.io/badge/R-4.0+-blue.svg)](https://www.r-project.org/)

## Overview

**MeLSI** (Metric Learning for Statistical Inference) is a novel machine learning method for microbiome data analysis that learns optimal distance metrics to improve statistical power in detecting group differences. Unlike traditional distance metrics (Bray-Curtis, Euclidean, Jaccard), MeLSI adapts to the specific characteristics of your dataset to maximize separation between groups.

If you use the MeLSI software, please cite our work:

Nathan Bresette. MeLSI: Metric Learning for Statistical Inference in Microbiome Analysis. 2025. https://github.com/NathanBresette/MeLSI

## Support

Check out the MeLSI examples in the `vignettes/` folder for an overview of analysis options and example runs. Users should start with the basic usage examples and then refer to the comprehensive testing suite as necessary.

If you have questions, please direct them to the [MeLSI Issues](https://github.com/NathanBresette/MeLSI/issues) page.

## Contents

- [Introduction](#introduction)
- [Support](#support)
- [Contents](#contents)
- [Requirements](#requirements)
- [Installation](#installation)
- [Running MeLSI](#running-melsi)
- [Input data](#input-data)
- [Output files](#output-files)
- [Run a demo](#run-a-demo)
- [Options](#options)
- [Validation Results](#validation-results)
- [Troubleshooting](#troubleshooting)

## Requirements

MeLSI is an R package that can be run from R or sourced directly.

## Installation

### Install from GitHub

The latest version of MeLSI can be installed from GitHub:

```r
# Install required packages first
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install Bioconductor packages
BiocManager::install(c("phyloseq", "microbiome"))

# Install CRAN packages
install.packages(c("vegan", "ggplot2", "dplyr", "GUniFrac"))

# Source MeLSI functions directly
source("R/melsi_robust.R")
```

## Running MeLSI

MeLSI can be run as an R function. The method requires microbiome count data, group labels, and optional parameters for customization. MeLSI will return statistical results including F-statistics, p-values, and learned metric weights.

### Input data

MeLSI requires microbiome count data formatted as:

**Feature abundance data frame**
- Formatted with features (taxa) as columns and samples as rows
- Contains raw count data or relative abundances
- Can include zeros (they will be handled appropriately)
- This can be a filepath to a tab-delimited file

**Group labels**
- Vector of group assignments for each sample
- Should be factors with meaningful group names
- Must match the order of samples in the abundance data

### Output files

MeLSI generates several types of output:

**Statistical Results**
- F-statistic: Test statistic for group differences
- P-value: Permutation-based p-value
- Learned metric: Matrix of learned feature weights
- Null distribution: F-statistics from permuted data
- Diagnostics: Method performance metrics

**Validation Results** (when running comprehensive tests)
- Performance comparisons with standard methods
- Type I error and power analysis results
- Scalability and parameter sensitivity results
- Biological validation results

## Run a demo

Example analysis can be found in the `R/` folder. The following code demonstrates basic MeLSI usage:

```r
# Source MeLSI functions
source("R/melsi_robust.R")

# Generate synthetic microbiome data
test_data <- generate_test_data(n_samples = 60, n_taxa = 100, n_signal_taxa = 10)
X <- test_data$counts
y <- test_data$metadata$Group

# CLR transformation (recommended)
X_clr <- X
X_clr[X_clr == 0] <- 1e-10
X_clr <- log(X_clr)
X_clr <- X_clr - rowMeans(X_clr)

# Run MeLSI analysis
results <- run_melsi_permutation_test(
    X_clr, y, 
    n_perms = 99,    # Number of permutations
    B = 50,          # Number of weak learners
    m_frac = 0.7,    # Fraction of features per learner
    show_progress = TRUE
)

# Display results
cat("MeLSI Results:\n")
cat(sprintf("F-statistic: %.4f\n", results$F_observed))
cat(sprintf("P-value: %.4f\n", results$p_value))
cat(sprintf("Significant: %s\n", ifelse(results$p_value < 0.05, "Yes", "No")))
```

## Options

### Required parameters

- `X`: A matrix of feature abundances with samples as rows and features as columns
- `y`: A vector of group labels for each sample

### Analysis options

- `n_perms` (default 99): Number of permutations for p-value calculation
- `B` (default 50): Number of weak learners in the ensemble
- `m_frac` (default 0.7): Fraction of features to use in each weak learner
- `pre_filter` (default TRUE): Whether to apply pre-filtering to remove low-variance features
- `show_progress` (default TRUE): Whether to display progress information

### Advanced options

- `learning_rate` (default 0.1): Learning rate for gradient descent optimization
- `filter_threshold` (default 0.1): Threshold for pre-filtering feature selection
- `max_iterations` (default 50): Maximum iterations for weak learner optimization

## Validation Results

### Performance Comparison

MeLSI was tested against 6 standard methods on multiple datasets:

| Method | Average F-statistic | Detection Rate |
|--------|-------------------|----------------|
| **MeLSI** | **2.45** | **80%** |
| Bray-Curtis | 1.68 | 60% |
| Euclidean | 1.89 | 67% |
| Jaccard | 1.09 | 33% |
| Weighted UniFrac | 1.32 | 40% |
| Unweighted UniFrac | 1.34 | 40% |

### Statistical Validation

- ✅ **Type I Error**: Maintains ~5% error rate under null hypothesis
- ✅ **Power Analysis**: Superior detection across effect sizes (60-70% improvement)
- ✅ **Scalability**: Efficient up to 1000+ taxa and 500+ samples
- ✅ **Biological Relevance**: 70% overlap with traditional differential abundance

### Real Data Validation

**Atlas1006 Sex Comparison**:
- MeLSI F-statistic: 6.05 vs Euclidean: 4.73 (28% improvement)
- Successfully detected known biological patterns

**SoilRep Warming Study**:
- MeLSI successfully detected environmental warming effects
- Robust performance on challenging high-dimensional soil microbiome data

## Troubleshooting

**Question**: When I run MeLSI I see convergence warnings. Is this normal?
**Answer**: Some convergence warnings are normal, especially with small datasets. Check the diagnostics output for condition numbers and eigenvalue ranges.

**Question**: MeLSI is running slowly on my large dataset. How can I speed it up?
**Answer**: Reduce the number of permutations (`n_perms`), ensemble size (`B`), or enable pre-filtering (`pre_filter = TRUE`).

**Question**: How do I interpret the learned metric weights?
**Answer**: Higher weights indicate taxa that contribute more to group separation. The diagonal elements of the learned metric matrix represent feature importance.

**Question**: Should I use CLR transformation?
**Answer**: Yes, CLR transformation is recommended for microbiome data as it handles compositionality and zeros appropriately.

## Citation

If you use MeLSI in your research, please cite:

```bibtex
@software{melsi2025,
  title={MeLSI: Metric Learning for Statistical Inference in Microbiome Analysis},
  author={Bresette, Nathan},
  year={2025},
  url={https://github.com/NathanBresette/MeLSI},
  note={Novel machine learning method for microbiome differential analysis}
}
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

*MeLSI: Advancing microbiome analysis through adaptive metric learning*
