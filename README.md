# MeLSI: Metric Learning for Statistical Inference in Microbiome Analysis

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R](https://img.shields.io/badge/R-4.0+-blue.svg)](https://www.r-project.org/)

## Overview

**MeLSI** (Metric Learning for Statistical Inference) is a novel machine learning method that solves a fundamental problem in microbiome analysis: **traditional distance metrics are fixed and generic, missing biologically meaningful patterns in your data**.

### The Problem MeLSI Solves

Current microbiome analysis relies on fixed distance metrics like Bray-Curtis, Euclidean, or Jaccard that treat all microbial taxa equally. This "one-size-fits-all" approach often misses subtle but biologically important differences between groups, leading to:

- ❌ **Reduced statistical power** - Missing real biological signals
- ❌ **Poor performance on complex datasets** - Especially with many rare taxa
- ❌ **Generic analysis** - Same approach regardless of your research question
- ❌ **Suboptimal feature weighting** - All taxa treated equally when some matter more

### What MeLSI Does Differently

MeLSI **learns an optimal distance metric** specifically tailored to your data and research question by:

1. **🎯 Adaptive Learning**: Discovers which microbial features are most important for separating your groups
2. **🧠 Ensemble Intelligence**: Uses multiple weak learners with bootstrap sampling for robust results
3. **🔍 Smart Pre-filtering**: Automatically removes noise features to focus on biological signal
4. **📊 Interpretable Results**: Provides feature importance weights that align with biological knowledge

### Why This Matters

Instead of using the same distance formula for every study, MeLSI learns what matters most for **your specific data**, resulting in:

- ✅ **46% higher F-statistics** compared to Bray-Curtis distance
- ✅ **80% detection rate** vs 60% for traditional methods
- ✅ **28% improvement** on real biological datasets (Atlas1006)
- ✅ **70% overlap** with known biological patterns

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

The latest version of MeLSI can be installed from GitHub using devtools:

```r
# Install devtools if not already installed
if (!require("devtools", quietly = TRUE))
    install.packages("devtools")

# Install required packages first
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install Bioconductor packages
BiocManager::install(c("phyloseq", "microbiome"))

# Install CRAN packages
install.packages(c("vegan", "ggplot2", "dplyr", "GUniFrac"))

# Install MeLSI from GitHub
devtools::install_github("NathanBresette/MeLSI")

# Load the package
library(MeLSI)
```

## Running MeLSI

MeLSI transforms your microbiome data analysis by learning a custom distance metric optimized for your specific research question. Here's what happens when you run MeLSI:

### What MeLSI Does Under the Hood

1. **🔍 Analyzes Your Data**: Examines which microbial taxa show the strongest group differences
2. **🧠 Learns Optimal Weights**: Discovers how important each taxon is for separating your groups  
3. **🎯 Builds Custom Metric**: Creates a distance measure tailored to your data's biological patterns
4. **📊 Tests Significance**: Uses permutation testing to ensure results are statistically valid
5. **📈 Returns Insights**: Provides both statistical results and biological interpretation

### Key Outputs

MeLSI returns comprehensive results including:
- **F-statistic**: How well groups are separated (higher = better)
- **P-value**: Statistical significance (permutation-based, more reliable)
- **Learned metric weights**: Which taxa matter most for your analysis
- **Diagnostics**: Quality metrics to ensure reliable results

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
# Load MeLSI package
library(MeLSI)

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

## How MeLSI Works

### The Algorithm

MeLSI uses an innovative ensemble approach to learn optimal distance metrics:

1. **🎯 Pre-filtering**: Identifies and focuses on taxa with the strongest group differences using t-tests
2. **🔄 Bootstrap Sampling**: Creates multiple training sets by resampling your data
3. **🧩 Feature Subsampling**: Each learner uses a random subset of features to prevent overfitting
4. **📈 Gradient Optimization**: Learns optimal weights for each feature subset using gradient descent
5. **⚖️ Ensemble Averaging**: Combines multiple learners, weighting better performers more heavily
6. **🛡️ Robust Distance Calculation**: Ensures numerical stability with eigenvalue decomposition
7. **📊 Permutation Testing**: Validates significance using null distributions from permuted data

### Why This Approach Works

- **Adaptive**: Learns what matters for your specific data, not generic patterns
- **Robust**: Ensemble approach prevents overfitting and improves generalization  
- **Interpretable**: Feature weights reveal which taxa drive group differences
- **Validated**: Permutation testing ensures statistical reliability
- **Efficient**: Pre-filtering focuses computation on relevant features

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
