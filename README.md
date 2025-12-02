# MeLSI: Metric Learning for Statistical Inference in Microbiome Analysis

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R](https://img.shields.io/badge/R-4.0+-blue.svg)](https://www.r-project.org/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17714848.svg)](https://doi.org/10.5281/zenodo.17714848)
[![Paper](https://img.shields.io/badge/Paper-mSystems-blue)](https://journals.asm.org/journal/msystems)

## Overview

**MeLSI** (Metric Learning for Statistical Inference) is a novel machine learning framework for microbiome beta diversity analysis that solves a fundamental problem: **traditional distance metrics are fixed and generic, missing biologically meaningful patterns in your data**. MeLSI learns adaptive distance metrics optimized for your specific dataset while maintaining rigorous statistical validity through permutation testing.

### The Problem

Current microbiome beta diversity analysis relies on fixed distance metrics like Bray-Curtis, Euclidean, or Jaccard that treat all microbial taxa equally. This "one-size-fits-all" approach often misses subtle but biologically important differences between groups.

### How MeLSI Works

MeLSI uses an innovative ensemble approach to learn optimal distance metrics for community composition analysis:

1. **Conservative Pre-filtering**: Selects taxa with highest variance using variance-based filtering, keeping 70% of features to maintain statistical power while reducing noise
2. **Bootstrap Sampling**: Creates multiple training sets by resampling your data
3. **Feature Subsampling**: Each learner uses a random subset of features (80%) to prevent overfitting
4. **Gradient Optimization**: Learns optimal weights for each feature subset using gradient descent
5. **Ensemble Averaging**: Combines 30 weak learners, weighting better performers more heavily
6. **Robust Distance Calculation**: Ensures numerical stability with eigenvalue decomposition
7. **Permutation Testing**: Validates significance using null distributions from permuted data (200 permutations for reliable p-values)

### Performance Results

Instead of using the same distance formula for every study, MeLSI learns what matters most for **your specific data**, resulting in:

- ✅ **Proper Type I Error Control** - No false positives on null data (p = 0.607 and 0.224)
- ✅ **Competitive Statistical Power** - Context-dependent performance with interpretable feature weights
- ✅ **Biological Interpretability** - Feature importance weights reveal which taxa drive group differences
- ✅ **Real-World Validation** - 9.1% improvement on Atlas1006, significant detection on DietSwap where traditional methods were marginal

## Requirements

### System Requirements
- **R version**: 4.0.0 or higher
- **Operating System**: Windows, macOS, or Linux

### Required R Packages
MeLSI depends on the following R packages (automatically installed):

**Core Dependencies:**
- `vegan` - Community ecology package for PERMANOVA calculations
- `ggplot2` - Graphics and visualization
- `dplyr` - Data manipulation
- `GUniFrac` - Phylogenetic distance calculations
- `stats` - Statistical functions (base R)
- `methods` - S4 methods (base R)

**Optional Dependencies:**
- `phyloseq` - Microbiome data structures (via BiocManager)
- `microbiome` - Microbiome analysis tools (via BiocManager)
- `devtools` - For GitHub installation
- `BiocManager` - For Bioconductor packages

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

### Key Outputs

MeLSI automatically adapts to your data:
- **2 groups**: Standard pairwise analysis with feature importance
- **3+ groups**: Comprehensive omnibus and pairwise analysis with multiple testing correction

Results include:
- **F-statistic**: How well groups are separated (higher = better)
- **P-value**: Statistical significance (permutation-based, more reliable)
- **Feature importance weights**: Which taxa matter most - automatically displayed as VIP charts
- **Directionality information**: Indicates which group has higher abundance for each important taxon (NEW!)
- **Log2 fold-change values**: Quantitative differences between groups (NEW!)
- **Top features**: Displayed in console with their learned weights and directionality
- **Multiple testing correction**: Automatic FDR correction for multi-group comparisons
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

## Quick Start (3 lines of code!)

The easiest way to use MeLSI:

```r
library(MeLSI)

# Your microbiome data: X = counts (samples × taxa), y = group labels
X_clr <- clr_transform(X)  # CLR transformation (recommended)
results <- melsi(X_clr, y)   # Run analysis - VIP plot auto-generated with directionality!
plot_pcoa(results, X_clr, y)  # PCoA plot

# That's it! Results include F-statistic, p-value, and feature importance
```

### Key improvements for ease of use:
- `clr_transform()` — One-line CLR transformation
- `plot_vip(results)` — Auto-includes directionality coloring (no extra parameters!)
- `plot_pcoa(results, X, y)` — Easy PCoA plotting from results
- Directionality information automatically calculated and displayed

## Detailed Examples

### Pairwise Analysis (2 groups)
```r
# Load MeLSI package
library(MeLSI)

# Generate synthetic microbiome data with realistic taxonomic names
set.seed(42)
n_samples <- 60
n_taxa <- 50

# Create realistic bacterial species names
realistic_taxa <- c(
    "Bacteroides_vulgatus", "Bacteroides_thetaiotaomicron", "Bacteroides_fragilis",
    "Faecalibacterium_prausnitzii", "Ruminococcus_bromii", "Ruminococcus_gnavus",
    "Bifidobacterium_longum", "Lactobacillus_acidophilus", "Akkermansia_muciniphila",
    "Prevotella_copri", "Parabacteroides_distasonis", "Alistipes_putredinis",
    "Roseburia_intestinalis", "Coprococcus_comes", "Blautia_wexlerae",
    "Collinsella_aerofaciens", "Eggerthella_lenta", "Fusobacterium_nucleatum",
    "Veillonella_parvula", "Streptococcus_salivarius", "Enterococcus_faecalis",
    "Staphylococcus_epidermidis", "Propionibacterium_acnes", "Corynebacterium_striatum",
    "Escherichia_coli", "Klebsiella_pneumoniae", "Enterobacter_cloacae",
    "Clostridium_innocuum", "Clostridium_leptum", "Methanobrevibacter_smithii",
    "Desulfovibrio_piger", "Succinivibrio_dextrinosolvens", "Micrococcus_luteus",
    "Sarcina_ventriculi", "Slackia_equolifaciens", "Barnesiella_intestinihominis",
    "Odoribacter_splanchnicus", "Porphyromonas_gingivalis", "Bacteroides_ovatus",
    "Bacteroides_uniformis", "Bacteroides_caccae", "Bifidobacterium_adolescentis",
    "Bifidobacterium_bifidum", "Lactobacillus_casei", "Lactobacillus_plantarum",
    "Lactobacillus_rhamnosus", "Lactobacillus_reuteri", "Prevotella_oralis",
    "Dorea_longicatena", "Ruminococcus_bromii"
)

# Create base microbiome data (log-normal distribution)
X <- matrix(rlnorm(n_samples * n_taxa, meanlog = 2, sdlog = 1), 
            nrow = n_samples, ncol = n_taxa)

# IMPORTANT: Set realistic column names
colnames(X) <- realistic_taxa[1:n_taxa]

# Create group labels
y <- c(rep("Control", n_samples/2), rep("Treatment", n_samples/2))

# Add signal to first 10 taxa in Treatment group (simulate differential abundance)
X[31:60, 1:10] <- X[31:60, 1:10] * 1.5

# CLR transformation (recommended for microbiome data) - NOW EASIER!
X_clr <- clr_transform(X)

# Run MeLSI analysis - VIP plot with directionality is auto-generated!
results <- melsi(X_clr, y, n_perms = 200, B = 30, show_progress = TRUE)

# Display results summary
cat("MeLSI Results:\n")
cat(sprintf("F-statistic: %.4f\n", results$F_observed))
cat(sprintf("P-value: %.4f\n", results$p_value))
cat(sprintf("Significant: %s\n", ifelse(results$p_value < 0.05, "Yes", "No")))

# Create PCoA plot - NOW SUPER EASY!
plot_pcoa(results, X_clr, y)

# Re-plot VIP with more features (directionality is automatic!)
plot_vip(results, top_n = 20)

# Access directionality information
cat("\nTop taxa with directionality:\n")
if (!is.null(results$directionality)) {
    dir_df <- data.frame(
        Taxon = names(results$directionality),
        Higher_in = results$directionality,
        Weight = results$feature_weights[names(results$directionality)],
        Log2FC = results$log2_fold_change
    )
    print(head(dir_df[order(-dir_df$Weight), ], 10))
}
```

### Multi-Group Analysis (3+ groups)
```r
# Create 4-group data with realistic taxonomic names
set.seed(42)
n_samples_per_group <- 30
n_taxa <- 50

# Use same realistic taxa as pairwise example
realistic_taxa <- c(
    "Bacteroides_vulgatus", "Bacteroides_thetaiotaomicron", "Bacteroides_fragilis",
    "Faecalibacterium_prausnitzii", "Ruminococcus_bromii", "Ruminococcus_gnavus",
    "Bifidobacterium_longum", "Lactobacillus_acidophilus", "Akkermansia_muciniphila",
    "Prevotella_copri", "Parabacteroides_distasonis", "Alistipes_putredinis",
    "Roseburia_intestinalis", "Coprococcus_comes", "Blautia_wexlerae",
    "Collinsella_aerofaciens", "Eggerthella_lenta", "Fusobacterium_nucleatum",
    "Veillonella_parvula", "Streptococcus_salivarius", "Enterococcus_faecalis",
    "Staphylococcus_epidermidis", "Propionibacterium_acnes", "Corynebacterium_striatum",
    "Escherichia_coli", "Klebsiella_pneumoniae", "Enterobacter_cloacae",
    "Clostridium_innocuum", "Clostridium_leptum", "Methanobrevibacter_smithii",
    "Desulfovibrio_piger", "Succinivibrio_dextrinosolvens", "Micrococcus_luteus",
    "Sarcina_ventriculi", "Slackia_equolifaciens", "Barnesiella_intestinihominis",
    "Odoribacter_splanchnicus", "Porphyromonas_gingivalis", "Bacteroides_ovatus",
    "Bacteroides_uniformis", "Bacteroides_caccae", "Bifidobacterium_adolescentis",
    "Bifidobacterium_bifidum", "Lactobacillus_casei", "Lactobacillus_plantarum",
    "Lactobacillus_rhamnosus", "Lactobacillus_reuteri", "Prevotella_oralis",
    "Dorea_longicatena", "Ruminococcus_bromii"
)

X <- matrix(rlnorm(n_samples_per_group * 4 * n_taxa, meanlog = 2, sdlog = 1), 
            nrow = n_samples_per_group * 4, ncol = n_taxa)

# IMPORTANT: Set realistic column names
colnames(X) <- realistic_taxa[1:n_taxa]

y <- rep(c("Control", "Mild", "Moderate", "Severe"), each = n_samples_per_group)

# Add progressive signal to simulate disease severity
X[31:60, 1:10] <- X[31:60, 1:10] * 1.4    # Mild
X[61:90, 1:10] <- X[61:90, 1:10] * 1.7    # Moderate  
X[91:120, 1:10] <- X[91:120, 1:10] * 2.0  # Severe

# CLR transformation
X_clr <- log(X + 1)
X_clr <- X_clr - rowMeans(X_clr)

# IMPORTANT: Preserve column names after transformation
colnames(X_clr) <- colnames(X)

# Run MeLSI analysis (automatically detects 4 groups and runs both omnibus and pairwise)
results <- melsi(X_clr, y, n_perms = 200, B = 30, show_progress = TRUE)

# Access omnibus results
cat("Omnibus F-statistic:", results$omnibus$F_observed, "\n")
cat("Omnibus p-value:", results$omnibus$p_value, "\n")

# Show global feature importance
cat("\nTop 5 globally important taxa:\n")
global_top <- head(sort(results$omnibus$feature_weights, decreasing = TRUE), 5)
for (i in 1:length(global_top)) {
    cat(sprintf("  %d. %-30s (weight: %.4f)\n", i, names(global_top)[i], global_top[i]))
}

# Access pairwise results
cat("\nSignificant pairwise comparisons:", nrow(results$pairwise$significant_pairs), "\n")
print(results$pairwise$summary_table)
```

## Key Advantages

- **Adaptive**: Learns what matters for your specific data, not generic patterns
- **Robust**: Ensemble approach prevents overfitting and improves generalization  
- **Interpretable**: Feature weights reveal which taxa drive group differences - automatically visualized with VIP charts
- **Directional**: Automatically determines which group has higher abundance for each important taxon (NEW!)
- **Validated**: Permutation testing ensures statistical reliability and proper Type I error control
- **Efficient**: Pre-filtering focuses computation on relevant features
- **User-friendly**: Automatic feature importance visualization with directionality coloring shows exactly which taxa matter most and how they differ

## Options

### Required parameters

- `X`: A matrix of feature abundances with samples as rows and features as columns
- `y`: A vector of group labels for each sample

### Analysis options

- `analysis_type` (default "auto"): Type of analysis - "auto", "pairwise", "omnibus", or "both"
- `n_perms` (default 200): Number of permutations for p-value calculation
- `B` (default 30): Number of weak learners in the ensemble
- `m_frac` (default 0.8): Fraction of features to use in each weak learner
- `show_progress` (default TRUE): Whether to display progress information
- `plot_vip` (default TRUE): Whether to automatically display Variable Importance Plot
- `correction_method` (default "BH"): Multiple testing correction method for multi-group analysis

### Advanced options

- `learning_rate` (default 0.1): Learning rate for gradient descent optimization
- `filter_threshold` (default 0.1): Threshold for pre-filtering feature selection
- `max_iterations` (default 50): Maximum iterations for weak learner optimization

## Validation Results

*Note: All traditional methods used in comparisons are PERMANOVA tests with standard distance metrics. MeLSI uses 200 permutations; traditional methods use 999 permutations.*

### Type I Error Control

**Table 1. Type I Error Control (Proper Calibration)**

| Dataset | MeLSI F | MeLSI p | Best Traditional | Best Trad F | Best Trad p |
|---------|---------|---------|------------------|-------------|-------------|
| **Null Synthetic** | 1.307 | 0.607 | Euclidean | 0.964 | 0.638 |
| **Null Real Shuffled** | 1.737 | 0.224 | Bray-Curtis | 1.020 | 0.397 |

*✅ MeLSI shows proper Type I error control - no false positives on null data! All p-values are appropriately high, well above 0.05.*

### Performance Comparison

**Table 2. Method Comparison on Synthetic & Real Data**

| Dataset | MeLSI F | MeLSI p | Best Traditional | Best Trad F | Best Trad p |
|---------|---------|---------|------------------|-------------|-------------|
| **Synthetic Small (1.5×)** | 1.333 | 0.373 | Weighted UniFrac | 1.592 | 0.021* |
| **Synthetic Medium (2.0×)** | 1.605 | 0.030* | Bray-Curtis | 1.829 | 0.001* |
| **Synthetic Large (3.0×)** | 2.217 | 0.005* | Weighted UniFrac | 6.145 | 0.001* |
| **Atlas1006 (Real)** | 5.141 | 0.005* | Euclidean | 4.711 | 0.001* |
| **DietSwap (Real)** | 2.856 | 0.015* | Bray-Curtis | 2.153 | 0.058 |

(*p < 0.05)

**Key Findings:**
- **Atlas1006**: MeLSI achieved 9.1% improvement over best traditional method (Euclidean)
- **DietSwap**: MeLSI detected significant differences (p=0.015) where best traditional method was marginal (p=0.058)
- **Synthetic data**: MeLSI shows appropriate conservatism on weak signals, detects medium/large effects reliably
- **Context-dependent**: Specialized methods (Weighted UniFrac, Bray-Curtis) may outperform on large effect sizes with phylogenetic structure

### Scalability Analysis

**Table 3. Computational Performance Across Sample Sizes**

| n | p | MeLSI F | MeLSI Time (s) | Best Traditional | Best Trad F | Best Trad Time (s) |
|---|---|---------|----------------|------------------|-------------|-------------------|
| 20 | 200 | 1.222 | 185.4 | Bray-Curtis | 1.133 | 0.014 |
| 50 | 200 | 1.263 | 181.6 | Bray-Curtis | 1.222 | 0.029 |
| 100 | 200 | 1.510 | 238.2 | Bray-Curtis | 1.676 | 0.087 |
| 200 | 200 | 1.548 | 480.0 | Bray-Curtis | 2.254 | 0.311 |
| 500 | 200 | 2.424 | 2244.3 | Bray-Curtis | 4.319 | 2.324 |

**Scalability Notes:**
- MeLSI shows good small-sample properties (appropriate conservatism at n=20-50)
- F-statistics increase monotonically with sample size as expected
- Computation time scales approximately O(n²) but remains practical for typical microbiome studies
- Pre-filtering provides computational benefits on high-dimensional datasets

### Parameter Sensitivity

**Table 4. Robust Hyperparameter Performance**

| Parameter | Value | F-statistic | p-value | Time (s) |
|-----------|-------|-------------|---------|----------|
| **Ensemble Size (B)** |
| | 10 | 1.438 | 0.179 | 98.7 |
| | 20 | 1.467 | 0.109 | 160.8 |
| | 30 | 1.478 | 0.090 | 235.0 |
| | 50 | 1.465 | 0.119 | 389.9 |
| | 100 | 1.462 | 0.100 | 768.1 |
| **Feature Fraction (m_frac)** |
| | 0.5 | 1.492 | 0.139 | 187.2 |
| | 0.7 | 1.459 | 0.109 | 213.5 |
| | 0.8 | 1.442 | 0.134 | 240.7 |
| | 0.9 | 1.422 | 0.124 | 262.2 |
| | 1.0 | 1.427 | 0.124 | 283.7 |

**Robustness Notes:**
- F-statistics remain stable across wide range of hyperparameters
- Default settings (B=30, m_frac=0.8) provide good balance of performance and computational efficiency
- Results validate that MeLSI is not overly sensitive to parameter choices

### Pre-filtering Benefits

**Table 5. Impact of Conservative Pre-filtering**

| Dataset | Effect | Features | With Filter F | Without Filter F | F Change | Time Saved |
|---------|--------|----------|---------------|------------------|----------|------------|
| Test 1 | Small | 500 | 1.278 | 1.284 | -0.5% | 5.8% |
| Test 2 | Medium | 200 | 1.432 | 1.416 | +1.7% | 4.1% |
| Test 3 | Large | 100 | 1.224 | 1.267 | -4.3% | 1.2% |

**Pre-filtering Notes:**
- Modest computational benefits (1-6% time reduction)
- Minimal impact on statistical power (changes < 5%)
- Benefits are context-dependent; most valuable for extremely high-dimensional datasets
- Conservative threshold (70% retention) prioritizes maintaining signal

## Reproducibility

All results in the research paper are fully reproducible. The `reproducibility_scripts/` directory contains standalone R scripts for generating each table:

- `table1_type1_error.R` - Type I error control analysis
- `table2_power_analysis.R` - Method comparison on synthetic and real datasets
- `table2_dietswap.R` - DietSwap-specific analysis
- `table3_scalability.R` - Scalability analysis across sample sizes
- `table4_parameter_sensitivity.R` - Parameter sensitivity analysis
- `table5_prefiltering.R` - Pre-filtering benefit analysis
- `run_all_tables.R` - Master script to run all analyses

Each script is self-contained, uses `set.seed(42)` for reproducibility, and generates CSV outputs. See `reproducibility_scripts/README.md` for detailed instructions.

## Troubleshooting

**Question**: When I run MeLSI I see convergence warnings. Is this normal?
**Answer**: Some convergence warnings are normal, especially with small datasets. Check the diagnostics output for condition numbers and eigenvalue ranges.

**Question**: MeLSI is running slowly on my large dataset. How can I speed it up?
**Answer**: Reduce the number of permutations (`n_perms`), ensemble size (`B`), or enable pre-filtering (`pre_filter = TRUE`).

**Question**: How do I interpret the learned metric weights?
**Answer**: Higher weights indicate taxa that contribute more to group separation. MeLSI automatically displays the top 5 most important features and generates VIP charts showing feature importance. For multi-group analysis, you get both global feature importance (omnibus) and comparison-specific importance (pairwise). Access all weights via `results$feature_weights` or `results$omnibus$feature_weights` and `results$pairwise$pairwise_results`.

**Question**: Should I use CLR transformation?
**Answer**: Yes, CLR transformation is recommended for microbiome data as it handles compositionality appropriately. MeLSI is designed to work with CLR-transformed data.

**Question**: Why does MeLSI sometimes give higher p-values than traditional methods?
**Answer**: This is actually a feature, not a bug! MeLSI is appropriately conservative and avoids false positives on borderline effects, while still detecting strong signals reliably. The interpretable feature weights provide biological insight even when p-values are similar to traditional methods.

**Question**: When should I choose MeLSI over traditional methods?
**Answer**: Prefer MeLSI for PERMANOVA workflows that need calibrated p-values plus interpretable taxa weights (e.g., diet interventions or subtle host phenotype comparisons). Traditional methods may be preferable when you need the fastest possible computation or when working with very large effect sizes and phylogenetic structure.

## Support

Check out the MeLSI examples in the `vignettes/` folder for an overview of analysis options and example runs. Users should start with the basic usage examples and then refer to the comprehensive testing suite as necessary.

If you have questions, please direct them to the [MeLSI Issues](https://github.com/NathanBresette/MeLSI/issues) page.

## Citation

If you use MeLSI in your research, please cite:

**Manuscript (In Preparation for mSystems):**
```bibtex
@article{bresette2025melsi,
  title={MeLSI: Metric Learning for Statistical Inference in Microbiome Community Composition Analysis},
  author={Bresette, Nathan and Ericsson, Aaron C. and Woods, Carter and Lin, Ai-Ling},
  journal={mSystems},
  year={2025},
  note={In preparation},
  doi={[To be added upon publication]}
}
```

**Software:**
```bibtex
@software{melsi_software,
  title={MeLSI: Metric Learning for Statistical Inference in Microbiome Analysis},
  author={Bresette, Nathan and Ericsson, Aaron C. and Woods, Carter and Lin, Ai-Ling},
  year={2025},
  url={https://github.com/NathanBresette/MeLSI},
  doi={10.5281/zenodo.17714848},
  note={Novel machine learning framework for microbiome beta diversity analysis}
}
```

*Note: DOIs will be updated upon publication and Zenodo archival*

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## Manuscript

The MeLSI methodology and comprehensive validation are described in detail in our manuscript currently in preparation for mSystems (American Society for Microbiology). The manuscript includes:

- Complete mathematical framework and algorithmic details
- Comprehensive validation across synthetic and real datasets
- Type I error control and statistical power analysis
- Scalability, parameter sensitivity, and computational performance evaluation
- Biological interpretability case studies with Atlas1006 and DietSwap datasets

All analyses presented in the manuscript are fully reproducible using the scripts provided in the `reproducibility_scripts/` directory.

---

*MeLSI: Advancing microbiome beta diversity analysis through adaptive metric learning, interpretable feature importance, and rigorous statistical inference*
