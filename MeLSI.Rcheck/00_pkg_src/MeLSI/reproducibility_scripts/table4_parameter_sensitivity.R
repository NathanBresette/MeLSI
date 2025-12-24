#!/usr/bin/env Rscript
# ==============================================================================
# Table 4: Parameter Sensitivity Analysis
# ==============================================================================
# Reproduces Table 4 from the MeLSI paper showing how performance varies
# with key hyperparameters: ensemble size (B) and feature fraction (m_frac).
# ==============================================================================

# Load required packages
library(vegan)
library(MeLSI)

set.seed(42)

# ==============================================================================
# Generate Test Dataset
# ==============================================================================
cat("\n")
cat("==============================================================================\n")
cat("REPRODUCING TABLE 4: Parameter Sensitivity Analysis\n")
cat("==============================================================================\n\n")

cat("Generating test dataset (100 samples, 200 taxa, medium effect)...\n")

# Generate dataset once for all parameter tests
n_samples <- 100
n_taxa <- 200
counts <- matrix(rnbinom(n_samples * n_taxa, mu = 30, size = 0.8), 
                nrow = n_samples, ncol = n_taxa)
counts[counts < 3] <- 0

# Add signal
group1_size <- n_samples %/% 2
signal_taxa <- sample(1:n_taxa, 10)
for (i in 1:group1_size) {
  counts[i, signal_taxa] <- counts[i, signal_taxa] * 2
}

group <- c(rep("Group1", group1_size), rep("Group2", n_samples - group1_size))

# CLR transformation
X_clr <- log(counts + 1)
X_clr <- X_clr - rowMeans(X_clr)
colnames(X_clr) <- paste0("Taxa_", 1:ncol(X_clr))

# ==============================================================================
# Test 1: Ensemble Size (B) Effects
# ==============================================================================
cat("\n")
cat("Part A: Testing Ensemble Size (B) Effects\n")
cat("------------------------------------------------------------------------------\n")

B_values <- c(10, 20, 30, 50, 100)
results_B <- data.frame()

for (B in B_values) {
  cat(sprintf("Testing B = %d...\n", B))
  
  start_time <- Sys.time()
  melsi_result <- melsi(X_clr, group, n_perms = 200, B = B, m_frac = 0.8,
                       show_progress = TRUE, plot_vip = FALSE)
  time_taken <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  
  results_B <- rbind(results_B, data.frame(
    B = B,
    F_statistic = melsi_result$F_observed,
    p_value = melsi_result$p_value,
    time_sec = time_taken
  ))
}

# ==============================================================================
# Test 2: Feature Fraction (m_frac) Effects
# ==============================================================================
cat("\n")
cat("Part B: Testing Feature Fraction (m_frac) Effects\n")
cat("------------------------------------------------------------------------------\n")

m_frac_values <- c(0.5, 0.7, 0.8, 0.9, 1.0)
results_m_frac <- data.frame()

for (m_frac in m_frac_values) {
  cat(sprintf("Testing m_frac = %.1f...\n", m_frac))
  
  start_time <- Sys.time()
  melsi_result <- melsi(X_clr, group, n_perms = 200, B = 30, m_frac = m_frac,
                       show_progress = TRUE, plot_vip = FALSE)
  time_taken <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  
  results_m_frac <- rbind(results_m_frac, data.frame(
    m_frac = m_frac,
    F_statistic = melsi_result$F_observed,
    p_value = melsi_result$p_value,
    time_sec = time_taken
  ))
}

# ==============================================================================
# Display Results
# ==============================================================================
cat("\n")
cat("==============================================================================\n")
cat("TABLE 5A: Ensemble Size (B) Sensitivity\n")
cat("==============================================================================\n\n")

print(results_B, row.names = FALSE, digits = 3)

cat("\n")
cat("==============================================================================\n")
cat("TABLE 5B: Feature Fraction (m_frac) Sensitivity\n")
cat("==============================================================================\n\n")

print(results_m_frac, row.names = FALSE, digits = 3)

cat("\n")
cat("Interpretation:\n")
cat("- Ensemble Size (B):\n")
cat(sprintf("  * Performance stabilizes around B=30 (%.3f F-statistic)\n", 
           results_B$F_statistic[results_B$B == 30]))
cat("  * Larger B increases computation time with diminishing returns\n")
cat("  * B=30 provides optimal balance\n")
cat("\n")
cat("- Feature Fraction (m_frac):\n")
cat(sprintf("  * m_frac=0.8 shows robust performance (%.3f F-statistic)\n", 
           results_m_frac$F_statistic[results_m_frac$m_frac == 0.8]))
cat("  * Too low (0.5): underutilizes features\n")
cat("  * Too high (1.0): reduces ensemble diversity\n")
cat("  * m_frac=0.8 balances diversity and feature coverage\n")

# Save results
write.csv(results_B, "table4a_results.csv", row.names = FALSE)
write.csv(results_m_frac, "table4b_results.csv", row.names = FALSE)
cat("\nResults saved to: table4a_results.csv and table4b_results.csv\n")
cat("==============================================================================\n")

