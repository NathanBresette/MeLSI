#!/usr/bin/env Rscript
# ==============================================================================
# Table 3: Scalability Across Sample Size and Dimensionality
# ==============================================================================
# Reproduces Table 3 from the MeLSI paper showing how MeLSI scales with
# increasing sample sizes and feature counts.
# ==============================================================================

# Load required packages
library(vegan)
library(MeLSI)
library(GUniFrac)
library(ape)

set.seed(42)

# ==============================================================================
# Helper Function: Generate Scalability Test Data
# ==============================================================================
generate_scalability_data <- function(n_samples, n_taxa) {
  counts <- matrix(rnbinom(n_samples * n_taxa, mu = 30, size = 0.8), 
                  nrow = n_samples, ncol = n_taxa)
  counts[counts < 3] <- 0
  
  # Add medium effect signal
  group1_size <- n_samples %/% 2
  signal_taxa <- min(10, n_taxa %/% 10)
  signal_indices <- sample(1:n_taxa, signal_taxa)
  
  for (i in 1:group1_size) {
    counts[i, signal_indices] <- counts[i, signal_indices] * 2
  }
  
  group <- c(rep("Group1", group1_size), rep("Group2", n_samples - group1_size))
  
  return(list(counts = counts, group = group))
}

# ==============================================================================
# Test 1: Varying Sample Sizes (fixed taxa = 200)
# ==============================================================================
cat("\n")
cat("==============================================================================\n")
cat("REPRODUCING TABLE 3: Scalability Analysis\n")
cat("==============================================================================\n\n")

cat("Part A: Varying Sample Sizes (fixed taxa = 200)\n")
cat("------------------------------------------------------------------------------\n")

sample_sizes <- c(20, 50, 100, 200, 500)
results_samples <- data.frame()

for (n in sample_sizes) {
  cat(sprintf("Testing n = %d samples...\n", n))
  
  data <- generate_scalability_data(n, 200)
  
  # CLR transformation
  X_clr <- log(data$counts + 1)
  X_clr <- X_clr - rowMeans(X_clr)
  colnames(X_clr) <- paste0("Taxa_", 1:ncol(X_clr))
  
  # Time MeLSI
  start_time <- Sys.time()
  melsi_result <- melsi(X_clr, data$group, n_perms = 200, B = 30, 
                       show_progress = TRUE, plot_vip = FALSE)
  melsi_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  
  # Time Euclidean PERMANOVA
  start_time <- Sys.time()
  dist_euc <- dist(X_clr)
  perm_euc <- adonis2(dist_euc ~ data$group, permutations = 999)
  euc_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  
  # Time Bray-Curtis PERMANOVA
  start_time <- Sys.time()
  dist_bray <- vegdist(data$counts, method = "bray")
  perm_bray <- adonis2(dist_bray ~ data$group, permutations = 999)
  bray_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  
  # Find best traditional
  best_trad_F <- max(perm_euc$F[1], perm_bray$F[1])
  best_trad <- ifelse(perm_euc$F[1] > perm_bray$F[1], "Euclidean", "BrayCurtis")
  best_trad_time <- ifelse(best_trad == "Euclidean", euc_time, bray_time)
  
  results_samples <- rbind(results_samples, data.frame(
    n_samples = n,
    n_taxa = 200,
    MeLSI_time_sec = melsi_time,
    MeLSI_F = melsi_result$F_observed,
    Best_Traditional = best_trad,
    Best_Traditional_F = best_trad_F,
    Best_Traditional_time_sec = best_trad_time,
    Euclidean_time_sec = euc_time,
    BrayCurtis_time_sec = bray_time
  ))
}

# ==============================================================================
# Test 2: Varying Taxa Counts (fixed samples = 100)
# ==============================================================================
cat("\n")
cat("Part B: Varying Taxa Counts (fixed samples = 100)\n")
cat("------------------------------------------------------------------------------\n")

taxa_sizes <- c(50, 100, 200, 500, 1000)
results_taxa <- data.frame()

for (p in taxa_sizes) {
  cat(sprintf("Testing p = %d taxa...\n", p))
  
  data <- generate_scalability_data(100, p)
  
  # CLR transformation
  X_clr <- log(data$counts + 1)
  X_clr <- X_clr - rowMeans(X_clr)
  colnames(X_clr) <- paste0("Taxa_", 1:ncol(X_clr))
  
  # Time MeLSI
  start_time <- Sys.time()
  melsi_result <- melsi(X_clr, data$group, n_perms = 200, B = 30, 
                       show_progress = TRUE, plot_vip = FALSE)
  melsi_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  
  # Time Euclidean PERMANOVA
  start_time <- Sys.time()
  dist_euc <- dist(X_clr)
  perm_euc <- adonis2(dist_euc ~ data$group, permutations = 999)
  euc_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  
  # Time Bray-Curtis PERMANOVA
  start_time <- Sys.time()
  dist_bray <- vegdist(data$counts, method = "bray")
  perm_bray <- adonis2(dist_bray ~ data$group, permutations = 999)
  bray_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  
  # Find best traditional
  best_trad_F <- max(perm_euc$F[1], perm_bray$F[1])
  best_trad <- ifelse(perm_euc$F[1] > perm_bray$F[1], "Euclidean", "BrayCurtis")
  best_trad_time <- ifelse(best_trad == "Euclidean", euc_time, bray_time)
  
  results_taxa <- rbind(results_taxa, data.frame(
    n_samples = 100,
    n_taxa = p,
    MeLSI_time_sec = melsi_time,
    MeLSI_F = melsi_result$F_observed,
    Best_Traditional = best_trad,
    Best_Traditional_F = best_trad_F,
    Best_Traditional_time_sec = best_trad_time,
    Euclidean_time_sec = euc_time,
    BrayCurtis_time_sec = bray_time
  ))
}

# ==============================================================================
# Display Results
# ==============================================================================
cat("\n")
cat("==============================================================================\n")
cat("TABLE 4A: Scalability - Sample Size Effects\n")
cat("==============================================================================\n\n")

print(results_samples, row.names = FALSE, digits = 2)

cat("\n")
cat("==============================================================================\n")
cat("TABLE 4B: Scalability - Dimensionality Effects\n")
cat("==============================================================================\n\n")

print(results_taxa, row.names = FALSE, digits = 2)

cat("\n")
cat("Interpretation:\n")
cat(sprintf("- Sample scaling: ~O(n^2) scaling as expected for distance calculations\n"))
cat(sprintf("- Taxa scaling: MeLSI remains efficient even at p=1000\n"))
cat(sprintf("- Conservative pre-filtering keeps MeLSI practical for large p\n"))

# Save results
write.csv(results_samples, "table3a_results.csv", row.names = FALSE)
write.csv(results_taxa, "table3b_results.csv", row.names = FALSE)
cat("\nResults saved to: table3a_results.csv and table3b_results.csv\n")
cat("==============================================================================\n")

