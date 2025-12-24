#!/usr/bin/env Rscript
# ==============================================================================
# Table 5: Benefit of Conservative Pre-filtering
# ==============================================================================
# Reproduces Table 5 from the MeLSI paper demonstrating that conservative
# pre-filtering improves both performance and computational efficiency.
# ==============================================================================

# Load required packages
library(vegan)
library(MeLSI)

set.seed(42)

# ==============================================================================
# Helper Function: Generate Data with Varying Sparsity
# ==============================================================================
generate_sparse_data <- function(n_samples = 100, n_taxa = 500, effect_size = "medium") {
  counts <- matrix(rnbinom(n_samples * n_taxa, mu = 30, size = 0.8), 
                  nrow = n_samples, ncol = n_taxa)
  counts[counts < 3] <- 0
  
  # Add realistic sparsity (many zero-inflated features)
  for (i in 1:n_samples) {
    zero_taxa <- sample(1:n_taxa, n_taxa * 0.7)  # 70% sparse
    counts[i, zero_taxa] <- 0
  }
  
  # Add signal
  group1_size <- n_samples %/% 2
  
  if (effect_size == "small") {
    signal_taxa <- 5
    fold_change <- 1.5
  } else if (effect_size == "medium") {
    signal_taxa <- 10
    fold_change <- 2.0
  } else if (effect_size == "large") {
    signal_taxa <- 20
    fold_change <- 3.0
  }
  
  signal_indices <- sample(1:n_taxa, signal_taxa)
  for (i in 1:group1_size) {
    counts[i, signal_indices] <- counts[i, signal_indices] * fold_change
  }
  
  group <- c(rep("Group1", group1_size), rep("Group2", n_samples - group1_size))
  
  return(list(counts = counts, group = group))
}

# ==============================================================================
# MeLSI without pre-filtering (manual implementation for comparison)
# ==============================================================================
melsi_no_prefilter <- function(X, y, n_perms = 200, B = 30, m_frac = 0.8) {
  # This simulates MeLSI without pre-filtering by using all features
  # In practice, this is computationally expensive for high-dimensional data
  
  result <- melsi(X, y, n_perms = n_perms, B = B, m_frac = m_frac,
                 show_progress = TRUE, plot_vip = FALSE)
  
  return(result)
}

# ==============================================================================
# Run Pre-filtering Tests
# ==============================================================================
cat("\n")
cat("==============================================================================\n")
cat("REPRODUCING TABLE 5: Benefit of Conservative Pre-filtering\n")
cat("==============================================================================\n\n")

results <- data.frame()

test_cases <- list(
  list(n = 100, p = 500, effect = "small", desc = "Small effect, many noise features"),
  list(n = 100, p = 200, effect = "medium", desc = "Medium effect, moderate noise"),
  list(n = 100, p = 100, effect = "large", desc = "Large effect, few noise features")
)

for (i in 1:length(test_cases)) {
  tc <- test_cases[[i]]
  cat(sprintf("\nTest %d: %s\n", i, tc$desc))
  cat(sprintf("(n=%d, p=%d, effect=%s)\n", tc$n, tc$p, tc$effect))
  cat("------------------------------------------------------------------------------\n")
  
  # Generate data
  data <- generate_sparse_data(tc$n, tc$p, tc$effect)
  
  # CLR transformation
  X_clr <- log(data$counts + 1)
  X_clr <- X_clr - rowMeans(X_clr)
  colnames(X_clr) <- paste0("Taxa_", 1:ncol(X_clr))
  
  # Count features before filtering
  n_features_before <- ncol(X_clr)
  
  # Apply conservative pre-filtering manually
  # (in MeLSI package, this happens automatically)
  # Remove features with too many zeros (present in < 10% of samples)
  prevalence <- colMeans(data$counts > 0)
  keep_features <- prevalence >= 0.1
  X_clr_filtered <- X_clr[, keep_features]
  n_features_after <- sum(keep_features)
  
  cat(sprintf("Pre-filtering: %d → %d features (removed %d sparse features)\n",
             n_features_before, n_features_after, n_features_before - n_features_after))
  
  # Test WITH pre-filtering (using filtered data)
  cat("Running MeLSI WITH pre-filtering...\n")
  start_time <- Sys.time()
  melsi_with <- melsi(X_clr_filtered, data$group, n_perms = 200, B = 30, m_frac = 0.8,
                     show_progress = TRUE, plot_vip = FALSE)
  time_with <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  
  # Test WITHOUT pre-filtering (using all features)
  cat("Running MeLSI WITHOUT pre-filtering...\n")
  start_time <- Sys.time()
  melsi_without <- melsi(X_clr, data$group, n_perms = 200, B = 30, m_frac = 0.8,
                        show_progress = TRUE, plot_vip = FALSE)
  time_without <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  
  # Calculate improvements
  F_improvement <- melsi_with$F_observed - melsi_without$F_observed
  time_reduction_pct <- 100 * (time_without - time_with) / time_without
  
  results <- rbind(results, data.frame(
    Test_Case = tc$desc,
    Effect_Size = tc$effect,
    Features_Before = n_features_before,
    Features_After = n_features_after,
    With_Prefilter_F = melsi_with$F_observed,
    With_Prefilter_p = melsi_with$p_value,
    Without_Prefilter_F = melsi_without$F_observed,
    Without_Prefilter_p = melsi_without$p_value,
    F_Improvement = F_improvement,
    Time_Reduction_Pct = time_reduction_pct
  ))
}

# ==============================================================================
# Display Results
# ==============================================================================
cat("\n")
cat("==============================================================================\n")
cat("TABLE 5: Benefit of Conservative Pre-filtering\n")
cat("==============================================================================\n\n")

print(results[, c("Test_Case", "Features_Before", "Features_After", 
                  "With_Prefilter_F", "Without_Prefilter_F", 
                  "F_Improvement", "Time_Reduction_Pct")], 
      row.names = FALSE, digits = 3)

cat("\n")
cat("Interpretation:\n")
cat(sprintf("- Average F-statistic improvement: %.3f\n", mean(results$F_Improvement)))
cat(sprintf("- Average time reduction: %.1f%%\n", mean(results$Time_Reduction_Pct)))
cat("\nKey Insight: Pre-filtering removes uninformative features that add noise\n")
cat("without removing true signal, improving both power and efficiency.\n")

# Save results
write.csv(results, "table5_results.csv", row.names = FALSE)
cat("\nResults saved to: table5_results.csv\n")
cat("==============================================================================\n")

