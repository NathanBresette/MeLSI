## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE)

## ----install, eval=FALSE------------------------------------------------------
# # Install devtools if not already installed
# if (!require("devtools", quietly = TRUE))
#     install.packages("devtools")
# 
# # Install required packages
# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# 
# BiocManager::install(c("phyloseq", "microbiome"))
# install.packages(c("vegan", "ggplot2"))
# 
# # Install MeLSI from GitHub
# devtools::install_github("NathanBresette/MeLSI")

## ----basic_setup--------------------------------------------------------------
# Load required packages
library(MeLSI)
library(vegan)
library(ggplot2)

# Generate synthetic microbiome data for demonstration
set.seed(42)
test_data <- generate_test_data(n_samples = 60, n_taxa = 100, n_signal_taxa = 10)
X <- test_data$counts
y <- test_data$metadata$Group

# Display data dimensions
cat("Data dimensions:", nrow(X), "samples x", ncol(X), "taxa\n")
cat("Groups:", paste(unique(y), collapse = ", "), "\n")

## ----preprocessing------------------------------------------------------------
# CLR transformation (recommended for microbiome data)
# Use the built-in helper function for easy CLR transformation
X_clr <- clr_transform(X)

# Display preprocessing results
cat("CLR transformation completed\n")
cat("Data range:", round(range(X_clr), 3), "\n")

## ----run_melsi----------------------------------------------------------------
# Run MeLSI with default parameters
# Note: For real analyses, use n_perms = 200+ for reliable p-values
# Using 99 here for faster vignette execution
cat("Running MeLSI analysis...\n")
melsi_results <- melsi(
    X_clr, y, 
    n_perms = 99,    # Number of permutations (use 200+ for real analysis)
    B = 30,          # Number of weak learners
    m_frac = 0.8,    # Fraction of features per learner
    show_progress = TRUE
)

# Display results
cat("\nMeLSI Results:\n")
cat(sprintf("F-statistic: %.4f\n", melsi_results$F_observed))
cat(sprintf("P-value: %.4f\n", melsi_results$p_value))
cat(sprintf("Significant: %s\n", ifelse(melsi_results$p_value < 0.05, "Yes", "No")))

# VIP plot is automatically generated, but you can also create it manually:
plot_vip(melsi_results, top_n = 15)

# PCoA plot using the learned MeLSI distance
plot_pcoa(melsi_results, X_clr, y)

## ----comparison---------------------------------------------------------------
# Compare with Euclidean distance
dist_euclidean <- dist(X_clr)
adonis_euclidean <- adonis2(dist_euclidean ~ y, permutations = 99)

# Compare with Bray-Curtis distance
dist_bray <- vegdist(X, method = "bray")
adonis_bray <- adonis2(dist_bray ~ y, permutations = 99)

# Display comparison results
cat("\nMethod Comparison:\n")
cat(sprintf("MeLSI:        F = %.4f, p = %.4f\n", 
           melsi_results$F_observed, melsi_results$p_value))
cat(sprintf("Euclidean:    F = %.4f, p = %.4f\n", 
           adonis_euclidean$F[1], adonis_euclidean$`Pr(>F)`[1]))
cat(sprintf("Bray-Curtis: F = %.4f, p = %.4f\n", 
           adonis_bray$F[1], adonis_bray$`Pr(>F)`[1]))

# Calculate improvement
improvement <- (melsi_results$F_observed - adonis_euclidean$F[1]) / adonis_euclidean$F[1] * 100
cat(sprintf("\nMeLSI improvement over Euclidean: %.1f%%\n", improvement))

## ----vip_plot-----------------------------------------------------------------
# Plot VIP with directionality (default)
plot_vip(melsi_results, top_n = 15, title = "Top 15 Features by Importance")

# Plot VIP without directionality coloring
plot_vip(melsi_results, top_n = 15, directionality = FALSE)

# Extract and examine feature weights programmatically
feature_weights <- melsi_results$feature_weights
top_features <- head(sort(feature_weights, decreasing = TRUE), 10)

cat("\nTop 10 Most Weighted Features:\n")
for (i in 1:length(top_features)) {
    cat(sprintf("%2d. %s (weight: %.4f)\n", 
               i, names(top_features)[i], top_features[i]))
}

## ----pcoa_plot----------------------------------------------------------------
# Create PCoA plot using the learned MeLSI distance matrix
plot_pcoa(melsi_results, X_clr, y, title = "PCoA: MeLSI Distance")

## ----advanced_usage-----------------------------------------------------------
# Run MeLSI with custom parameters
cat("Running MeLSI with custom parameters...\n")
custom_results <- melsi(
    X_clr, y, 
    n_perms = 199,   # More permutations for higher precision (recommended: 200+)
    B = 100,         # Larger ensemble for more stable results
    m_frac = 0.8,    # Fraction of features per learner
    show_progress = TRUE
)

cat(sprintf("Custom MeLSI F-statistic: %.4f (p = %.4f)\n", 
           custom_results$F_observed, custom_results$p_value))

## ----session_info-------------------------------------------------------------
sessionInfo()

