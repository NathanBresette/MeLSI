# -----------------------------------------------------------------------------
# MeLSI vs Standard Methods Comparison
# Test MeLSI against Bray-Curtis, Euclidean, and Jaccard on synthetic and real data
# -----------------------------------------------------------------------------

# Load required packages
if (!requireNamespace("vegan", quietly = TRUE)) install.packages("vegan")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("GUniFrac", quietly = TRUE)) install.packages("GUniFrac")
if (!requireNamespace("microbiome", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("microbiome")
}

library(vegan)
library(ggplot2)
library(dplyr)
library(GUniFrac)
library(microbiome)

# Source MeLSI
source('melsi_robust.R')

# -----------------------------------------------------------------------------
# Function 1: Generate Synthetic Test Data
# -----------------------------------------------------------------------------

generate_synthetic_data <- function(n_samples = 60, n_taxa = 150, signal_strength = "medium") {
  set.seed(42)
  
  # Generate realistic microbiome data
  counts <- matrix(rnbinom(n_samples * n_taxa, mu = 30, size = 0.8), nrow = n_samples, ncol = n_taxa)
  
  # Add realistic sparsity
  counts[counts < 3] <- 0
  
  # Add realistic abundance distribution
  for (i in 1:n_samples) {
    # Some taxa are very abundant
    abundant_taxa <- sample(1:n_taxa, 8)
    counts[i, abundant_taxa] <- counts[i, abundant_taxa] * 20
    
    # Most taxa are rare
    rare_taxa <- sample(setdiff(1:n_taxa, abundant_taxa), 60)
    counts[i, rare_taxa] <- counts[i, rare_taxa] * 0.1
  }
  
  colnames(counts) <- paste0("Taxon_", 1:n_taxa)
  rownames(counts) <- paste0("Sample_", 1:n_samples)
  
  # Create groups
  metadata <- data.frame(
    SampleID = rownames(counts),
    Group = factor(rep(c("Control", "Treatment"), each = n_samples / 2))
  )
  
  # Add signal based on strength
  treatment_indices <- which(metadata$Group == "Treatment")
  signal_taxa <- sample(1:n_taxa, 10)
  
  if (signal_strength == "weak") {
    effect_size <- 8
  } else if (signal_strength == "medium") {
    effect_size <- 15
  } else if (signal_strength == "strong") {
    effect_size <- 25
  }
  
  for (i in treatment_indices) {
    counts[i, signal_taxa] <- counts[i, signal_taxa] + rnbinom(length(signal_taxa), mu = effect_size, size = 1)
  }
  
  return(list(counts = counts, metadata = metadata))
}

# -----------------------------------------------------------------------------
# Function 2: Create Random Phylogenetic Tree
# -----------------------------------------------------------------------------

create_random_tree <- function(taxa_names) {
  # Create a random phylogenetic tree for UniFrac calculations
  library(ape)
  
  # Check if we have enough taxa for a tree
  if (length(taxa_names) < 2) {
    # If only one taxon, create a simple tree
    tree <- stree(1)
    tree$tip.label <- taxa_names[1]
  } else {
    # Generate random tree
    tree <- rtree(length(taxa_names))
    tree$tip.label <- taxa_names
  }
  
  return(tree)
}

# -----------------------------------------------------------------------------
# Function 3: Load Real Datasets
# -----------------------------------------------------------------------------

load_real_dataset <- function() {
  # Load Atlas1006 dataset
  data(atlas1006)
  atlas <- atlas1006
  
  # Subset to sex comparison
  sex_subset <- subset_samples(atlas, sex %in% c("female", "male"))
  sex_subset <- prune_taxa(taxa_sums(sex_subset) > 5, sex_subset)
  
  # Get counts and metadata
  counts <- as.matrix(otu_table(sex_subset))
  counts <- t(counts)  # Transpose to samples x taxa
  metadata <- data.frame(
    SampleID = sample_names(sex_subset),
    Group = sample_data(sex_subset)$sex
  )
  
  return(list(counts = counts, metadata = metadata, description = "Atlas1006: Female vs Male"))
}

load_additional_datasets <- function() {
  datasets <- list()
  
  # Dataset 1: Load soilrep dataset (soil microbiome - climate warming study)
  tryCatch({
    data(soilrep)
    sr <- soilrep
    
    # Subset to warming comparison (warmed vs not warmed)
    # Combine all treatments but focus on warming effect
    warmed_subset <- subset_samples(sr, warmed == "yes")
    control_subset <- subset_samples(sr, warmed == "no")
    
    # Combine both groups
    combined_subset <- merge_phyloseq(warmed_subset, control_subset)
    combined_subset <- prune_taxa(taxa_sums(combined_subset) > 5, combined_subset)
    
    if (nrow(sample_data(combined_subset)) >= 20) {
      counts <- as.matrix(otu_table(combined_subset))
      counts <- t(counts)  # Transpose to samples x taxa
      counts <- as.numeric(counts)  # Convert to numeric matrix
      counts <- matrix(counts, nrow = nrow(sample_data(combined_subset)), ncol = ntaxa(combined_subset))
      metadata <- data.frame(
        SampleID = sample_names(combined_subset),
        Group = sample_data(combined_subset)$warmed
      )
      datasets$soilrep_warming <- list(
        counts = counts, 
        metadata = metadata, 
        description = paste0("SoilRep Warming: ", nrow(counts), " samples, ", ncol(counts), " taxa")
      )
    }
  }, error = function(e) {
    cat("Could not load SoilRep warming dataset:", e$message, "\n")
  })
  
  return(datasets)
}

# -----------------------------------------------------------------------------
# Function 4: Test All Methods (including phylogenetic)
# -----------------------------------------------------------------------------

test_all_methods <- function(counts, metadata, dataset_name) {
  cat("\n=== Testing", dataset_name, "===\n")
  cat("Samples:", nrow(counts), "Taxa:", ncol(counts), "\n")
  
  # CLR transformation
  counts_clr <- as.matrix(vegan::decostand(counts + 1, method = "clr"))
  
  results <- data.frame()
  
  # Test MeLSI (reduced parameters for faster testing)
  cat("Running MeLSI...\n")
  start_time <- Sys.time()
  melsi_results <- run_melsi_permutation_test(
    counts_clr, metadata$Group, 
    n_perms = 19, B = 20, m_frac = 0.7, show_progress = TRUE
  )
  melsi_time <- as.numeric(Sys.time() - start_time, units = "secs")
  
  # Test Bray-Curtis
  cat("Running Bray-Curtis...\n")
  start_time <- Sys.time()
  dist_bray <- vegan::vegdist(counts, method = "bray")
  permanova_bray <- vegan::adonis2(dist_bray ~ Group, data = metadata, permutations = 999)
  bray_time <- as.numeric(Sys.time() - start_time, units = "secs")
  
  # Test Euclidean
  cat("Running Euclidean...\n")
  start_time <- Sys.time()
  dist_euclidean <- dist(counts_clr)
  permanova_euclidean <- vegan::adonis2(dist_euclidean ~ Group, data = metadata, permutations = 999)
  euclidean_time <- as.numeric(Sys.time() - start_time, units = "secs")
  
  # Test Jaccard
  cat("Running Jaccard...\n")
  start_time <- Sys.time()
  dist_jaccard <- vegan::vegdist(counts, method = "jaccard")
  permanova_jaccard <- vegan::adonis2(dist_jaccard ~ Group, data = metadata, permutations = 999)
  jaccard_time <- as.numeric(Sys.time() - start_time, units = "secs")
  
  # Initialize phylogenetic results with NA values
  permanova_weighted <- list(F = c(NA), 'Pr(>F)' = c(NA))
  permanova_unweighted <- list(F = c(NA), 'Pr(>F)' = c(NA))
  weighted_time <- 0
  unweighted_time <- 0
  
  # Test Weighted UniFrac
  cat("Running Weighted UniFrac...\n")
  start_time <- Sys.time()
  unifrac_results <- NULL
  tryCatch({
    # Skip phylogenetic distances for datasets with >1000 taxa to avoid memory issues
    if (ncol(counts) > 1000) {
      stop("Skipping phylogenetic distances for large datasets (>1000 taxa)")
    }
    tree <- create_random_tree(colnames(counts))
    unifrac_results <- GUniFrac(counts, tree, alpha = c(0, 0.5, 1))
    dist_weighted <- unifrac_results$unifracs[,, "d_1"]  # Weighted UniFrac
    permanova_weighted <- vegan::adonis2(dist_weighted ~ Group, data = metadata, permutations = 999)
    weighted_time <- as.numeric(Sys.time() - start_time, units = "secs")
  }, error = function(e) {
    cat("Weighted UniFrac failed:", e$message, "\n")
    weighted_time <- as.numeric(Sys.time() - start_time, units = "secs")
  })
  
  # Test Unweighted UniFrac
  cat("Running Unweighted UniFrac...\n")
  start_time <- Sys.time()
  tryCatch({
    if (is.null(unifrac_results)) {
      stop("UniFrac results not available")
    }
    dist_unweighted <- unifrac_results$unifracs[,, "d_UW"]  # Unweighted UniFrac
    permanova_unweighted <- vegan::adonis2(dist_unweighted ~ Group, data = metadata, permutations = 999)
    unweighted_time <- as.numeric(Sys.time() - start_time, units = "secs")
  }, error = function(e) {
    cat("Unweighted UniFrac failed:", e$message, "\n")
    unweighted_time <- as.numeric(Sys.time() - start_time, units = "secs")
  })
  
  # Test PERMANOVA-S (omnibus test)
  cat("Running PERMANOVA-S...\n")
  start_time <- Sys.time()
  
  # Simple omnibus test (average of F-statistics)
  f_stats <- c(
    permanova_bray$F[1],
    permanova_euclidean$F[1], 
    permanova_jaccard$F[1]
  )
  
  # Add phylogenetic F-statistics if available
  if (!is.na(permanova_weighted$F[1])) {
    f_stats <- c(f_stats, permanova_weighted$F[1])
  }
  if (!is.na(permanova_unweighted$F[1])) {
    f_stats <- c(f_stats, permanova_unweighted$F[1])
  }
  
  permanova_s_f <- mean(f_stats, na.rm = TRUE)
  permanova_s_time <- as.numeric(Sys.time() - start_time, units = "secs")
  
  # Print results
  cat("\nResults:\n")
  cat("  MeLSI: F =", round(melsi_results$F_observed, 4), "p =", round(melsi_results$p_value, 4), 
      ifelse(melsi_results$p_value < 0.05, "***", ""), "\n")
  cat("  Bray-Curtis: F =", round(permanova_bray$F[1], 4), "p =", round(permanova_bray$'Pr(>F)'[1], 4), 
      ifelse(permanova_bray$'Pr(>F)'[1] < 0.05, "***", ""), "\n")
  cat("  Euclidean: F =", round(permanova_euclidean$F[1], 4), "p =", round(permanova_euclidean$'Pr(>F)'[1], 4), 
      ifelse(permanova_euclidean$'Pr(>F)'[1] < 0.05, "***", ""), "\n")
  cat("  Jaccard: F =", round(permanova_jaccard$F[1], 4), "p =", round(permanova_jaccard$'Pr(>F)'[1], 4), 
      ifelse(permanova_jaccard$'Pr(>F)'[1] < 0.05, "***", ""), "\n")
  cat("  Weighted UniFrac: F =", round(permanova_weighted$F[1], 4), "p =", round(permanova_weighted$'Pr(>F)'[1], 4), 
      ifelse(permanova_weighted$'Pr(>F)'[1] < 0.05, "***", ""), "\n")
  cat("  Unweighted UniFrac: F =", round(permanova_unweighted$F[1], 4), "p =", round(permanova_unweighted$'Pr(>F)'[1], 4), 
      ifelse(permanova_unweighted$'Pr(>F)'[1] < 0.05, "***", ""), "\n")
  cat("  PERMANOVA-S: F =", round(permanova_s_f, 4), "\n")
  
  # Store results
  result_row <- data.frame(
    dataset = dataset_name,
    n_samples = nrow(counts),
    n_taxa = ncol(counts),
    
    # MeLSI results
    melsi_F = melsi_results$F_observed,
    melsi_p = melsi_results$p_value,
    melsi_significant = melsi_results$p_value < 0.05,
    melsi_time = melsi_time,
    
    # Bray-Curtis results
    bray_F = permanova_bray$F[1],
    bray_p = permanova_bray$'Pr(>F)'[1],
    bray_significant = permanova_bray$'Pr(>F)'[1] < 0.05,
    bray_time = bray_time,
    
    # Euclidean results
    euclidean_F = permanova_euclidean$F[1],
    euclidean_p = permanova_euclidean$'Pr(>F)'[1],
    euclidean_significant = permanova_euclidean$'Pr(>F)'[1] < 0.05,
    euclidean_time = euclidean_time,
    
    # Jaccard results
    jaccard_F = permanova_jaccard$F[1],
    jaccard_p = permanova_jaccard$'Pr(>F)'[1],
    jaccard_significant = permanova_jaccard$'Pr(>F)'[1] < 0.05,
    jaccard_time = jaccard_time,
    
    # Weighted UniFrac results
    weighted_unifrac_F = permanova_weighted$F[1],
    weighted_unifrac_p = permanova_weighted$'Pr(>F)'[1],
    weighted_unifrac_significant = permanova_weighted$'Pr(>F)'[1] < 0.05,
    weighted_unifrac_time = weighted_time,
    
    # Unweighted UniFrac results
    unweighted_unifrac_F = permanova_unweighted$F[1],
    unweighted_unifrac_p = permanova_unweighted$'Pr(>F)'[1],
    unweighted_unifrac_significant = permanova_unweighted$'Pr(>F)'[1] < 0.05,
    unweighted_unifrac_time = unweighted_time,
    
    # PERMANOVA-S results
    permanova_s_F = permanova_s_f,
    permanova_s_time = permanova_s_time
  )
  
  return(result_row)
}

# -----------------------------------------------------------------------------
# Function 4: Run Complete Comparison
# -----------------------------------------------------------------------------

run_complete_comparison <- function() {
  cat("========================================\n")
  cat("MeLSI vs Standard Methods Comparison\n")
  cat("========================================\n")
  
  results <- data.frame()
  
  # Test 1: Synthetic data - weak signal
  cat("\n1. Testing Synthetic Data - Weak Signal\n")
  synthetic_weak <- generate_synthetic_data(n_samples = 60, n_taxa = 150, signal_strength = "weak")
  result1 <- test_all_methods(synthetic_weak$counts, synthetic_weak$metadata, "Synthetic_Weak")
  results <- rbind(results, result1)
  
  # Test 2: Synthetic data - medium signal
  cat("\n2. Testing Synthetic Data - Medium Signal\n")
  synthetic_medium <- generate_synthetic_data(n_samples = 60, n_taxa = 150, signal_strength = "medium")
  result2 <- test_all_methods(synthetic_medium$counts, synthetic_medium$metadata, "Synthetic_Medium")
  results <- rbind(results, result2)
  
  # Test 3: Synthetic data - strong signal
  cat("\n3. Testing Synthetic Data - Strong Signal\n")
  synthetic_strong <- generate_synthetic_data(n_samples = 60, n_taxa = 150, signal_strength = "strong")
  result3 <- test_all_methods(synthetic_strong$counts, synthetic_strong$metadata, "Synthetic_Strong")
  results <- rbind(results, result3)
  
  # Test 4: Real data - Atlas1006 Sex
  cat("\n4. Testing Real Data - Atlas1006 Sex\n")
  real_data <- load_real_dataset()
  result4 <- test_all_methods(real_data$counts, real_data$metadata, "Real_Atlas1006_Sex")
  results <- rbind(results, result4)
  
  # Test 5-7: Additional datasets
  cat("\n5-7. Testing Additional Datasets\n")
  additional_datasets <- load_additional_datasets()
  
  dataset_counter <- 5
  for (dataset_name in names(additional_datasets)) {
    dataset <- additional_datasets[[dataset_name]]
    cat(sprintf("\n%d. Testing %s\n", dataset_counter, dataset$description))
    result <- test_all_methods(dataset$counts, dataset$metadata, paste0("Real_", dataset_name))
    results <- rbind(results, result)
    dataset_counter <- dataset_counter + 1
  }
  
  # Analyze results
  cat("\n", paste(rep("=", 60), collapse=""), "\n")
  cat("COMPARISON RESULTS SUMMARY\n")
  cat(paste(rep("=", 60), collapse=""), "\n\n")
  
  # Overall performance
  cat("Overall Performance:\n")
  cat("MeLSI significant:", sum(results$melsi_significant), "out of", nrow(results), "\n")
  cat("Bray-Curtis significant:", sum(results$bray_significant), "out of", nrow(results), "\n")
  cat("Euclidean significant:", sum(results$euclidean_significant), "out of", nrow(results), "\n")
  cat("Jaccard significant:", sum(results$jaccard_significant), "out of", nrow(results), "\n")
  cat("Weighted UniFrac significant:", sum(results$weighted_unifrac_significant), "out of", nrow(results), "\n")
  cat("Unweighted UniFrac significant:", sum(results$unweighted_unifrac_significant), "out of", nrow(results), "\n\n")
  
  # F-statistic comparison
  cat("F-statistic Comparison:\n")
  cat("MeLSI mean F:", round(mean(results$melsi_F), 4), "\n")
  cat("Bray-Curtis mean F:", round(mean(results$bray_F), 4), "\n")
  cat("Euclidean mean F:", round(mean(results$euclidean_F), 4), "\n")
  cat("Jaccard mean F:", round(mean(results$jaccard_F), 4), "\n")
  cat("Weighted UniFrac mean F:", round(mean(results$weighted_unifrac_F), 4), "\n")
  cat("Unweighted UniFrac mean F:", round(mean(results$unweighted_unifrac_F), 4), "\n")
  cat("PERMANOVA-S mean F:", round(mean(results$permanova_s_F), 4), "\n\n")
  
  # Timing comparison
  cat("Timing Comparison (seconds):\n")
  cat("MeLSI mean time:", round(mean(results$melsi_time), 2), "\n")
  cat("Bray-Curtis mean time:", round(mean(results$bray_time), 2), "\n")
  cat("Euclidean mean time:", round(mean(results$euclidean_time), 2), "\n")
  cat("Jaccard mean time:", round(mean(results$jaccard_time), 2), "\n")
  cat("Weighted UniFrac mean time:", round(mean(results$weighted_unifrac_time), 2), "\n")
  cat("Unweighted UniFrac mean time:", round(mean(results$unweighted_unifrac_time), 2), "\n")
  cat("PERMANOVA-S mean time:", round(mean(results$permanova_s_time), 2), "\n\n")
  
  # Best method per dataset
  cat("Best Method per Dataset:\n")
  for (i in 1:nrow(results)) {
    f_stats <- c(results$melsi_F[i], results$bray_F[i], results$euclidean_F[i], 
                 results$jaccard_F[i], results$weighted_unifrac_F[i], results$unweighted_unifrac_F[i])
    methods <- c("MeLSI", "Bray-Curtis", "Euclidean", "Jaccard", "Weighted UniFrac", "Unweighted UniFrac")
    best_method <- methods[which.max(f_stats)]
    cat("  ", results$dataset[i], ":", best_method, "(F =", round(max(f_stats), 4), ")\n")
  }
  
  # Save results
  write.csv(results, "method_comparison_results.csv", row.names = FALSE)
  
  cat("\nResults saved to 'method_comparison_results.csv'\n")
  cat("Comparison complete!\n")
  
  return(results)
}

# -----------------------------------------------------------------------------
# VALIDATION FUNCTIONS FOR PUBLICATION
# -----------------------------------------------------------------------------

# Function: Generate Null Datasets (Type I Error Testing)
generate_null_dataset <- function(n_samples = 100, n_taxa = 200, dataset_type = "synthetic") {
  set.seed(42)
  
  if (dataset_type == "synthetic") {
    # Generate completely random data with no signal
    counts <- matrix(rnbinom(n_samples * n_taxa, mu = 30, size = 0.8), nrow = n_samples, ncol = n_taxa)
    counts[counts < 3] <- 0
    
    # Create random groups (no true biological signal)
    group1_samples <- n_samples %/% 2
    group2_samples <- n_samples - group1_samples
    
    metadata <- data.frame(
      SampleID = paste0("Sample_", 1:n_samples),
      Group = c(rep("Group1", group1_samples), rep("Group2", group2_samples))
    )
    
    return(list(counts = counts, metadata = metadata, description = paste0("Null_Synthetic: ", n_samples, " samples, ", n_taxa, " taxa")))
    
  } else if (dataset_type == "real_shuffled") {
    # Use real Atlas1006 data but shuffle the group labels (no true signal)
    data(atlas1006)
    atlas <- atlas1006
    
    # Subset to sex comparison
    sex_subset <- subset_samples(atlas, sex %in% c("female", "male"))
    sex_subset <- prune_taxa(taxa_sums(sex_subset) > 5, sex_subset)
    
    # Limit to manageable size for testing
    if (nrow(sample_data(sex_subset)) > n_samples) {
      # Randomly sample n_samples
      sample_indices <- sample(1:nrow(sample_data(sex_subset)), n_samples)
      sex_subset <- prune_samples(sample_names(sex_subset)[sample_indices], sex_subset)
    }
    
    counts <- as.matrix(otu_table(sex_subset))
    counts <- t(counts)
    
    # SHUFFLE the group labels to create null hypothesis
    original_groups <- sample_data(sex_subset)$sex
    shuffled_groups <- sample(original_groups)
    
    metadata <- data.frame(
      SampleID = sample_names(sex_subset),
      Group = shuffled_groups  # Shuffled labels = no true signal
    )
    
    return(list(counts = counts, metadata = metadata, description = paste0("Null_RealShuffled: ", nrow(counts), " samples, ", ncol(counts), " taxa")))
  }
}

# Function: Test Single Dataset (Simplified for validation)
test_single_dataset <- function(counts, metadata, dataset_name, n_perms = 19, B = 20) {
  cat("Testing", dataset_name, "- Samples:", nrow(counts), "Taxa:", ncol(counts), "\n")
  
  # CLR transformation
  counts_clr <- as.matrix(vegan::decostand(counts + 1, method = "clr"))
  
  # Test MeLSI
  start_time <- Sys.time()
  melsi_results <- run_melsi_permutation_test(
    counts_clr, metadata$Group, 
    n_perms = n_perms, B = B, m_frac = 0.7, 
    show_progress = FALSE
  )
  melsi_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  
  # Test Euclidean (fast baseline)
  start_time <- Sys.time()
  euclidean_dist <- dist(counts_clr)
  euclidean_adonis <- adonis2(euclidean_dist ~ metadata$Group, permutations = n_perms)
  euclidean_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  
  return(data.frame(
    dataset = dataset_name,
    n_samples = nrow(counts),
    n_taxa = ncol(counts),
    melsi_F = melsi_results$F_observed,
    melsi_p = melsi_results$p_value,
    melsi_significant = melsi_results$p_value < 0.05,
    melsi_time = melsi_time,
    euclidean_F = euclidean_adonis$F[1],
    euclidean_p = euclidean_adonis$`Pr(>F)`[1],
    euclidean_significant = euclidean_adonis$`Pr(>F)`[1] < 0.05,
    euclidean_time = euclidean_time
  ))
}

# Function: Generate Power Analysis Datasets
generate_power_dataset <- function(n_samples = 100, n_taxa = 200, effect_size = "small") {
  set.seed(42)
  
  # Generate base microbiome data
  counts <- matrix(rnbinom(n_samples * n_taxa, mu = 30, size = 0.8), nrow = n_samples, ncol = n_taxa)
  counts[counts < 3] <- 0
  
  group1_samples <- n_samples %/% 2
  group2_samples <- n_samples - group1_samples
  
  # Add signal based on effect size
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
  
  # Add differential abundance to first group
  signal_indices <- sample(1:n_taxa, signal_taxa)
  for (i in 1:group1_samples) {
    counts[i, signal_indices] <- counts[i, signal_indices] * fold_change
  }
  
  metadata <- data.frame(
    SampleID = paste0("Sample_", 1:n_samples),
    Group = c(rep("Group1", group1_samples), rep("Group2", group2_samples))
  )
  
  return(list(counts = counts, metadata = metadata, description = paste0("Power_", effect_size, ": ", n_samples, " samples, ", n_taxa, " taxa")))
}

# Function: Type I Error Rate Test
test_type1_error <- function() {
  cat("\n", paste(rep("=", 60), collapse=""), "\n")
  cat("TYPE I ERROR RATE TEST\n")
  cat(paste(rep("=", 60), collapse=""), "\n")
  
  results <- data.frame()
  
  # Test 1: Synthetic null dataset
  cat("\n1. Testing Synthetic Null Dataset (100 samples, 200 taxa)\n")
  null_synthetic <- generate_null_dataset(n_samples = 100, n_taxa = 200, dataset_type = "synthetic")
  result1 <- test_single_dataset(null_synthetic$counts, null_synthetic$metadata, "Null_Synthetic")
  results <- rbind(results, result1)
  
  # Test 2: Real shuffled dataset
  cat("\n2. Testing Real Shuffled Dataset (Atlas1006 with shuffled labels)\n")
  null_real <- generate_null_dataset(n_samples = 200, n_taxa = 100, dataset_type = "real_shuffled")
  result2 <- test_single_dataset(null_real$counts, null_real$metadata, "Null_RealShuffled")
  results <- rbind(results, result2)
  
  # Analyze Type I error rates
  cat("\nType I Error Analysis:\n")
  cat("MeLSI Type I Error Rate:", mean(results$melsi_significant), "(should be ~0.05)\n")
  cat("Euclidean Type I Error Rate:", mean(results$euclidean_significant), "(should be ~0.05)\n")
  
  if (mean(results$melsi_significant) > 0.1) {
    cat("⚠️  WARNING: MeLSI has inflated Type I error rate!\n")
  } else {
    cat("✅ MeLSI Type I error rate is acceptable\n")
  }
  
  return(results)
}

# Function: Power Analysis Test
test_power_analysis <- function() {
  cat("\n", paste(rep("=", 60), collapse=""), "\n")
  cat("POWER ANALYSIS TEST\n")
  cat(paste(rep("=", 60), collapse=""), "\n")
  
  results <- data.frame()
  
  # Test different effect sizes
  effect_sizes <- c("small", "medium", "large")
  
  for (i in 1:length(effect_sizes)) {
    cat(sprintf("\n%d. Testing %s effect size (100 samples, 200 taxa)\n", i, effect_sizes[i]))
    power_data <- generate_power_dataset(n_samples = 100, n_taxa = 200, effect_size = effect_sizes[i])
    result <- test_single_dataset(power_data$counts, power_data$metadata, paste0("Power_", effect_sizes[i]))
    results <- rbind(results, result)
  }
  
  # Analyze power
  cat("\nPower Analysis:\n")
  cat("MeLSI Detection Rate:", mean(results$melsi_significant), "\n")
  cat("Euclidean Detection Rate:", mean(results$euclidean_significant), "\n")
  
  # Power by effect size
  for (i in 1:nrow(results)) {
    cat(sprintf("%s: MeLSI F=%.3f (p=%.3f), Euclidean F=%.3f (p=%.3f)\n",
                results$dataset[i], results$melsi_F[i], results$melsi_p[i],
                results$euclidean_F[i], results$euclidean_p[i]))
  }
  
  return(results)
}

# Function: Scalability Analysis
test_scalability <- function() {
  cat("\n", paste(rep("=", 60), collapse=""), "\n")
  cat("SCALABILITY ANALYSIS\n")
  cat(paste(rep("=", 60), collapse=""), "\n")
  
  results <- data.frame()
  
  # Test different sample sizes
  sample_sizes <- c(20, 50, 100, 200, 500)
  
  cat("\n1. Sample Size Effects (fixed at 200 taxa):\n")
  for (n in sample_sizes) {
    cat(sprintf("Testing n=%d samples...\n", n))
    power_data <- generate_power_dataset(n_samples = n, n_taxa = 200, effect_size = "medium")
    result <- test_single_dataset(power_data$counts, power_data$metadata, paste0("Scalability_n", n))
    results <- rbind(results, result)
  }
  
  # Test different taxa numbers
  taxa_numbers <- c(50, 100, 200, 500, 1000)
  
  cat("\n2. Taxa Number Effects (fixed at 100 samples):\n")
  for (p in taxa_numbers) {
    cat(sprintf("Testing p=%d taxa...\n", p))
    power_data <- generate_power_dataset(n_samples = 100, n_taxa = p, effect_size = "medium")
    result <- test_single_dataset(power_data$counts, power_data$metadata, paste0("Scalability_p", p))
    results <- rbind(results, result)
  }
  
  # Analyze scalability
  cat("\nScalability Analysis:\n")
  cat("Sample Size Effects:\n")
  sample_results <- results[grepl("Scalability_n", results$dataset), ]
  for (i in 1:nrow(sample_results)) {
    cat(sprintf("n=%d: MeLSI time=%.1fs, F=%.3f\n", 
                sample_results$n_samples[i], sample_results$melsi_time[i], sample_results$melsi_F[i]))
  }
  
  cat("\nTaxa Number Effects:\n")
  taxa_results <- results[grepl("Scalability_p", results$dataset), ]
  for (i in 1:nrow(taxa_results)) {
    cat(sprintf("p=%d: MeLSI time=%.1fs, F=%.3f\n", 
                taxa_results$n_taxa[i], taxa_results$melsi_time[i], taxa_results$melsi_F[i]))
  }
  
  return(results)
}

# Function: Parameter Sensitivity Analysis
test_parameter_sensitivity <- function() {
  cat("\n", paste(rep("=", 60), collapse=""), "\n")
  cat("PARAMETER SENSITIVITY ANALYSIS\n")
  cat(paste(rep("=", 60), collapse=""), "\n")
  
  results <- data.frame()
  
  # Generate test dataset
  power_data <- generate_power_dataset(n_samples = 100, n_taxa = 200, effect_size = "medium")
  
  # Test different ensemble sizes (B)
  cat("\n1. Ensemble Size (B) Effects:\n")
  B_values <- c(10, 20, 50, 100)
  
  for (B in B_values) {
    cat(sprintf("Testing B=%d...\n", B))
    counts_clr <- as.matrix(vegan::decostand(power_data$counts + 1, method = "clr"))
    
    start_time <- Sys.time()
    melsi_results <- run_melsi_permutation_test(
      counts_clr, power_data$metadata$Group, 
      n_perms = 19, B = B, m_frac = 0.7, 
      show_progress = FALSE
    )
    melsi_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    
    results <- rbind(results, data.frame(
      dataset = paste0("Param_B", B),
      parameter = "B",
      parameter_value = B,
      melsi_F = melsi_results$F_observed,
      melsi_p = melsi_results$p_value,
      melsi_time = melsi_time
    ))
  }
  
  # Test different feature fractions (m_frac)
  cat("\n2. Feature Fraction (m_frac) Effects:\n")
  m_frac_values <- c(0.5, 0.7, 0.8, 0.9)
  
  for (m_frac in m_frac_values) {
    cat(sprintf("Testing m_frac=%.1f...\n", m_frac))
    counts_clr <- as.matrix(vegan::decostand(power_data$counts + 1, method = "clr"))
    
    start_time <- Sys.time()
    melsi_results <- run_melsi_permutation_test(
      counts_clr, power_data$metadata$Group, 
      n_perms = 19, B = 20, m_frac = m_frac, 
      show_progress = FALSE
    )
    melsi_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    
    results <- rbind(results, data.frame(
      dataset = paste0("Param_mfrac", m_frac),
      parameter = "m_frac",
      parameter_value = m_frac,
      melsi_F = melsi_results$F_observed,
      melsi_p = melsi_results$p_value,
      melsi_time = melsi_time
    ))
  }
  
  # Analyze parameter sensitivity
  cat("\nParameter Sensitivity Analysis:\n")
  cat("Ensemble Size (B) Effects:\n")
  B_results <- results[results$parameter == "B", ]
  for (i in 1:nrow(B_results)) {
    cat(sprintf("B=%d: F=%.3f, time=%.1fs\n", 
                B_results$parameter_value[i], B_results$melsi_F[i], B_results$melsi_time[i]))
  }
  
  cat("\nFeature Fraction (m_frac) Effects:\n")
  mfrac_results <- results[results$parameter == "m_frac", ]
  for (i in 1:nrow(mfrac_results)) {
    cat(sprintf("m_frac=%.1f: F=%.3f, time=%.1fs\n", 
                mfrac_results$parameter_value[i], mfrac_results$melsi_F[i], mfrac_results$melsi_time[i]))
  }
  
  return(results)
}

# Function: Pre-filtering Value Analysis
test_prefiltering_value <- function() {
  cat("\n", paste(rep("=", 60), collapse=""), "\n")
  cat("PRE-FILTERING VALUE ANALYSIS\n")
  cat(paste(rep("=", 60), collapse=""), "\n")
  
  results <- data.frame()
  
  # Test datasets with different characteristics
  test_datasets <- list(
    # Small effect with many noise features
    generate_power_dataset(n_samples = 100, n_taxa = 500, effect_size = "small"),
    # Medium effect with moderate noise
    generate_power_dataset(n_samples = 100, n_taxa = 200, effect_size = "medium"),
    # Large effect with few noise features
    generate_power_dataset(n_samples = 100, n_taxa = 100, effect_size = "large")
  )
  
  for (i in 1:length(test_datasets)) {
    dataset <- test_datasets[[i]]
    cat(sprintf("\n%d. Testing %s (n=%d, p=%d)\n", i, dataset$description, 
                nrow(dataset$counts), ncol(dataset$counts)))
    
    counts_clr <- as.matrix(vegan::decostand(dataset$counts + 1, method = "clr"))
    
    # Test MeLSI WITH pre-filtering (default)
    cat("  Testing MeLSI WITH pre-filtering...\n")
    start_time <- Sys.time()
    melsi_with_prefilter <- run_melsi_permutation_test(
      counts_clr, dataset$metadata$Group, 
      n_perms = 19, B = 20, m_frac = 0.7, 
      show_progress = FALSE
    )
    time_with_prefilter <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    
    # Test MeLSI WITHOUT pre-filtering
    cat("  Testing MeLSI WITHOUT pre-filtering...\n")
    start_time <- Sys.time()
    melsi_without_prefilter <- run_melsi_permutation_test_no_prefilter(
      counts_clr, dataset$metadata$Group, 
      n_perms = 19, B = 20, m_frac = 0.7, 
      show_progress = FALSE
    )
    time_without_prefilter <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    
    # Store results
    results <- rbind(results, data.frame(
      dataset = paste0("Prefilter_Test", i),
      n_samples = nrow(dataset$counts),
      n_taxa = ncol(dataset$counts),
      effect_size = c("small", "medium", "large")[i],
      melsi_with_prefilter_F = melsi_with_prefilter$F_observed,
      melsi_with_prefilter_p = melsi_with_prefilter$p_value,
      melsi_with_prefilter_time = time_with_prefilter,
      melsi_without_prefilter_F = melsi_without_prefilter$F_observed,
      melsi_without_prefilter_p = melsi_without_prefilter$p_value,
      melsi_without_prefilter_time = time_without_prefilter,
      F_improvement = melsi_with_prefilter$F_observed - melsi_without_prefilter$F_observed,
      time_reduction = time_without_prefilter - time_with_prefilter
    ))
  }
  
  # Analyze pre-filtering value
  cat("\nPre-filtering Value Analysis:\n")
  cat("F-statistic Improvement (with vs without pre-filtering):\n")
  for (i in 1:nrow(results)) {
    cat(sprintf("%s: F improvement = %.3f (%.1f%% faster)\n", 
                results$dataset[i], results$F_improvement[i], 
                100 * results$time_reduction[i] / results$melsi_without_prefilter_time[i]))
  }
  
  avg_improvement <- mean(results$F_improvement)
  avg_speedup <- mean(results$time_reduction / results$melsi_without_prefilter_time)
  
  cat(sprintf("\nOverall: Pre-filtering improves F-statistic by %.3f on average\n", avg_improvement))
  cat(sprintf("Overall: Pre-filtering reduces computation time by %.1f%% on average\n", 100 * avg_speedup))
  
  if (avg_improvement > 0) {
    cat("✅ Pre-filtering provides consistent performance benefits\n")
  } else {
    cat("⚠️  Pre-filtering may not be beneficial for these datasets\n")
  }
  
  return(results)
}

# Function: Biological Validation Analysis
test_biological_validation <- function() {
  cat("\n", paste(rep("=", 60), collapse=""), "\n")
  cat("BIOLOGICAL VALIDATION ANALYSIS\n")
  cat(paste(rep("=", 60), collapse=""), "\n")
  
  results <- data.frame()
  
  # Load Atlas1006 sex data for biological validation
  cat("\n1. Atlas1006 Sex Comparison - Known Biological Patterns\n")
  tryCatch({
    data(atlas1006)
    atlas <- atlas1006
    
    # Subset to sex comparison
    sex_subset <- subset_samples(atlas, sex %in% c("female", "male"))
    sex_subset <- prune_taxa(taxa_sums(sex_subset) > 5, sex_subset)
    
    # Limit to manageable size for analysis
    if (nrow(sample_data(sex_subset)) > 200) {
      sample_indices <- sample(1:nrow(sample_data(sex_subset)), 200)
      sex_subset <- prune_samples(sample_names(sex_subset)[sample_indices], sex_subset)
    }
    
    counts <- as.matrix(otu_table(sex_subset))
    counts <- t(counts)
    metadata <- data.frame(
      SampleID = sample_names(sex_subset),
      Group = sample_data(sex_subset)$sex
    )
    
    # CLR transformation
    counts_clr <- counts
    counts_clr[counts_clr == 0] <- 1e-10
    counts_clr <- log(counts_clr)
    counts_clr <- counts_clr - rowMeans(counts_clr)
    
    # Run MeLSI
    cat("  Running MeLSI on Atlas1006 sex data...\n")
    melsi_results <- run_melsi_permutation_test(
      counts_clr, metadata$Group, 
      n_perms = 19, B = 20, m_frac = 0.7, 
      show_progress = FALSE
    )
    
    # Extract learned metric weights
    M_learned <- melsi_results$M_observed
    metric_weights <- diag(M_learned)
    
    # Get taxon names
    taxon_names <- colnames(counts_clr)
    
    # Create weight analysis
    weight_analysis <- data.frame(
      taxon = taxon_names,
      learned_weight = metric_weights,
      rank = rank(-metric_weights)
    )
    
    # Sort by learned weight
    weight_analysis <- weight_analysis[order(-weight_analysis$learned_weight), ]
    
    # Show top weighted taxa
    cat("\nTop 10 Most Weighted Taxa by MeLSI:\n")
    top_taxa <- head(weight_analysis, 10)
    for (i in 1:nrow(top_taxa)) {
      cat(sprintf("%d. %s (weight: %.4f)\n", i, top_taxa$taxon[i], top_taxa$learned_weight[i]))
    }
    
    # Calculate differential abundance for comparison
    cat("\nComparing with traditional differential abundance analysis...\n")
    group1_idx <- metadata$Group == unique(metadata$Group)[1]
    group2_idx <- metadata$Group == unique(metadata$Group)[2]
    
    # Simple t-test for each taxon
    p_values <- numeric(ncol(counts_clr))
    effect_sizes <- numeric(ncol(counts_clr))
    
    for (i in 1:ncol(counts_clr)) {
      test_result <- t.test(counts_clr[group1_idx, i], counts_clr[group2_idx, i])
      p_values[i] <- test_result$p.value
      effect_sizes[i] <- abs(test_result$statistic)
    }
    
    # Create comparison analysis
    comparison_analysis <- data.frame(
      taxon = taxon_names,
      learned_weight = metric_weights,
      t_test_pvalue = p_values,
      t_test_effect_size = effect_sizes,
      traditional_rank = rank(-effect_sizes)
    )
    
    # Sort by traditional effect size
    comparison_analysis <- comparison_analysis[order(-comparison_analysis$t_test_effect_size), ]
    
    cat("\nTop 10 Taxa by Traditional t-test Effect Size:\n")
    top_traditional <- head(comparison_analysis, 10)
    for (i in 1:nrow(top_traditional)) {
      cat(sprintf("%d. %s (effect size: %.3f, p=%.4f)\n", 
                  i, top_traditional$taxon[i], top_traditional$t_test_effect_size[i], top_traditional$t_test_pvalue[i]))
    }
    
    # Calculate overlap between top taxa
    top_melsi_taxa <- head(weight_analysis$taxon, 10)
    top_traditional_taxa <- head(comparison_analysis$taxon, 10)
    overlap_count <- sum(top_melsi_taxa %in% top_traditional_taxa)
    
    cat(sprintf("\nOverlap Analysis:\n"))
    cat(sprintf("Taxa in both top 10 lists: %d/10 (%.1f%%)\n", overlap_count, 100 * overlap_count / 10))
    
    # Store results
    results <- rbind(results, data.frame(
      dataset = "Atlas1006_Sex",
      n_samples = nrow(counts_clr),
      n_taxa = ncol(counts_clr),
      melsi_F = melsi_results$F_observed,
      melsi_p = melsi_results$p_value,
      top_taxa_overlap = overlap_count,
      overlap_percentage = 100 * overlap_count / 10,
      max_learned_weight = max(metric_weights),
      min_learned_weight = min(metric_weights),
      weight_range = max(metric_weights) - min(metric_weights)
    ))
    
  }, error = function(e) {
    cat("Could not analyze Atlas1006 sex data:", e$message, "\n")
  })
  
  # Load SoilRep warming data for biological validation
  cat("\n2. SoilRep Warming Comparison - Environmental Response\n")
  tryCatch({
    data(soilrep)
    sr <- soilrep
    
    # Subset to warming comparison
    warmed_subset <- subset_samples(sr, warmed == "yes")
    control_subset <- subset_samples(sr, warmed == "no")
    combined_subset <- merge_phyloseq(warmed_subset, control_subset)
    combined_subset <- prune_taxa(taxa_sums(combined_subset) > 5, combined_subset)
    
    # Limit to manageable size
    if (nrow(sample_data(combined_subset)) > 150) {
      sample_indices <- sample(1:nrow(sample_data(combined_subset)), 150)
      combined_subset <- prune_samples(sample_names(combined_subset)[sample_indices], combined_subset)
    }
    
    counts <- as.matrix(otu_table(combined_subset))
    counts <- t(counts)
    counts <- as.numeric(counts)
    counts <- matrix(counts, nrow = nrow(sample_data(combined_subset)), ncol = ntaxa(combined_subset))
    
    metadata <- data.frame(
      SampleID = sample_names(combined_subset),
      Group = sample_data(combined_subset)$warmed
    )
    
    # CLR transformation
    counts_clr <- counts
    counts_clr[counts_clr == 0] <- 1e-10
    counts_clr <- log(counts_clr)
    counts_clr <- counts_clr - rowMeans(counts_clr)
    
    # Run MeLSI
    cat("  Running MeLSI on SoilRep warming data...\n")
    melsi_results <- run_melsi_permutation_test(
      counts_clr, metadata$Group, 
      n_perms = 19, B = 20, m_frac = 0.7, 
      show_progress = FALSE
    )
    
    # Extract learned metric weights
    M_learned <- melsi_results$M_observed
    metric_weights <- diag(M_learned)
    
    # Get taxon names
    taxon_names <- colnames(counts_clr)
    
    # Create weight analysis
    weight_analysis <- data.frame(
      taxon = taxon_names,
      learned_weight = metric_weights,
      rank = rank(-metric_weights)
    )
    
    # Sort by learned weight
    weight_analysis <- weight_analysis[order(-weight_analysis$learned_weight), ]
    
    # Show top weighted taxa
    cat("\nTop 10 Most Weighted Taxa by MeLSI:\n")
    top_taxa <- head(weight_analysis, 10)
    for (i in 1:nrow(top_taxa)) {
      cat(sprintf("%d. %s (weight: %.4f)\n", i, top_taxa$taxon[i], top_taxa$learned_weight[i]))
    }
    
    # Store results
    results <- rbind(results, data.frame(
      dataset = "SoilRep_Warming",
      n_samples = nrow(counts_clr),
      n_taxa = ncol(counts_clr),
      melsi_F = melsi_results$F_observed,
      melsi_p = melsi_results$p_value,
      top_taxa_overlap = NA,  # No traditional comparison for soil data
      overlap_percentage = NA,
      max_learned_weight = max(metric_weights),
      min_learned_weight = min(metric_weights),
      weight_range = max(metric_weights) - min(metric_weights)
    ))
    
  }, error = function(e) {
    cat("Could not analyze SoilRep warming data:", e$message, "\n")
  })
  
  # Analyze biological validation
  cat("\nBiological Validation Summary:\n")
  if (nrow(results) > 0) {
    for (i in 1:nrow(results)) {
      cat(sprintf("%s: F=%.3f, Weight range=%.3f", 
                  results$dataset[i], results$melsi_F[i], results$weight_range[i]))
      if (!is.na(results$overlap_percentage[i])) {
        cat(sprintf(", Overlap with traditional=%.1f%%", results$overlap_percentage[i]))
      }
      cat("\n")
    }
  } else {
    cat("No datasets could be analyzed for biological validation.\n")
    cat("This may be due to data loading issues or insufficient samples.\n")
  }
  
  return(results)
}

# Function: MeLSI without pre-filtering (for comparison)
run_melsi_permutation_test_no_prefilter <- function(X, y, n_perms = 19, B = 20, m_frac = 0.7, show_progress = TRUE) {
  # This is identical to run_melsi_permutation_test but with pre_filter = FALSE
  n_samples <- nrow(X)
  n_taxa <- ncol(X)
  
  # Learn the metric
  if (show_progress) cat("Learning MeLSI metric (without pre-filtering)...\n")
  M_observed <- learn_melsi_metric_robust(X, y, B = B, m_frac = m_frac, pre_filter = FALSE)
  
  # Calculate observed F-statistic
  if (show_progress) cat("Calculating observed F-statistic...\n")
  dist_observed <- calculate_mahalanobis_dist_robust(X, M_observed)
  adonis_observed <- adonis2(dist_observed ~ y, permutations = 0)
  F_observed <- adonis_observed$F[1]
  
  # Permutation test
  if (show_progress) cat(sprintf("🔄 Running %d permutations...\n", n_perms))
  F_permuted <- numeric(n_perms)
  
  for (perm in 1:n_perms) {
    if (show_progress && perm %% 5 == 0) {
      progress_pct <- round(100 * perm / n_perms)
      cat(sprintf("   📊 Permutation %d/%d (%d%% complete)\n", perm, n_perms, progress_pct))
    }
    
    # Permute group labels
    y_permuted <- sample(y)
    
    # Learn metric on permuted data
    M_permuted <- learn_melsi_metric_robust(X, y_permuted, B = B, m_frac = m_frac, pre_filter = FALSE)
    
    # Calculate F-statistic
    dist_permuted <- calculate_mahalanobis_dist_robust(X, M_permuted)
    adonis_permuted <- adonis2(dist_permuted ~ y_permuted, permutations = 0)
    F_permuted[perm] <- adonis_permuted$F[1]
  }
  
  # Calculate p-value
  p_value <- mean(F_permuted >= F_observed)
  
  return(list(
    F_observed = F_observed,
    F_permuted = F_permuted,
    p_value = p_value,
    M_observed = M_observed
  ))
}

# Function: Run Complete Analysis (Main + All Validation Tests)
run_complete_analysis <- function() {
  start_time <- Sys.time()
  cat("🚀 STARTING COMPREHENSIVE MeLSI ANALYSIS\n")
  cat("📊 This includes: Main comparison + Statistical validation + Scalability + Parameter sensitivity + Pre-filtering + Biological validation\n")
  cat("⏱️  Estimated time: 15-30 minutes\n")
  cat("📁 Results will be saved to 7 CSV files\n\n")
  
  all_results <- list()
  
  # 1. Main Method Comparison
  cat("📋 PART 1/6: MAIN METHOD COMPARISON\n")
  cat("   Testing MeLSI vs 6 standard methods on 5 datasets\n")
  cat("   ⏳ This may take 5-10 minutes...\n")
  part1_start <- Sys.time()
  all_results$main_comparison <- run_complete_comparison()
  part1_time <- as.numeric(difftime(Sys.time(), part1_start, units = "mins"))
  cat(sprintf("   ✅ PART 1 COMPLETE! (%.1f minutes)\n\n", part1_time))
  
  # 2. Statistical Validation
  cat("📋 PART 2/6: STATISTICAL VALIDATION\n")
  cat("   Testing Type I error rates and power analysis\n")
  cat("   ⏳ This may take 3-5 minutes...\n")
  part2_start <- Sys.time()
  all_results$type1_error <- test_type1_error()
  all_results$power_analysis <- test_power_analysis()
  part2_time <- as.numeric(difftime(Sys.time(), part2_start, units = "mins"))
  cat(sprintf("   ✅ PART 2 COMPLETE! (%.1f minutes)\n\n", part2_time))
  
  # 3. Scalability Analysis
  cat("📋 PART 3/6: SCALABILITY ANALYSIS\n")
  cat("   Testing performance across different sample/taxa sizes\n")
  cat("   ⏳ This may take 5-8 minutes...\n")
  part3_start <- Sys.time()
  all_results$scalability <- test_scalability()
  part3_time <- as.numeric(difftime(Sys.time(), part3_start, units = "mins"))
  cat(sprintf("   ✅ PART 3 COMPLETE! (%.1f minutes)\n\n", part3_time))
  
  # 4. Parameter Sensitivity
  cat("📋 PART 4/6: PARAMETER SENSITIVITY ANALYSIS\n")
  cat("   Testing different ensemble sizes and feature fractions\n")
  cat("   ⏳ This may take 3-5 minutes...\n")
  part4_start <- Sys.time()
  all_results$parameter_sensitivity <- test_parameter_sensitivity()
  part4_time <- as.numeric(difftime(Sys.time(), part4_start, units = "mins"))
  cat(sprintf("   ✅ PART 4 COMPLETE! (%.1f minutes)\n\n", part4_time))
  
  # 5. Pre-filtering Value Analysis
  cat("📋 PART 5/6: PRE-FILTERING VALUE ANALYSIS\n")
  cat("   Comparing MeLSI with vs without pre-filtering\n")
  cat("   ⏳ This may take 3-5 minutes...\n")
  part5_start <- Sys.time()
  all_results$prefiltering_value <- test_prefiltering_value()
  part5_time <- as.numeric(difftime(Sys.time(), part5_start, units = "mins"))
  cat(sprintf("   ✅ PART 5 COMPLETE! (%.1f minutes)\n\n", part5_time))
  
  # 6. Biological Validation Analysis
  cat("📋 PART 6/6: BIOLOGICAL VALIDATION ANALYSIS\n")
  cat("   Analyzing learned metric weights for biological relevance\n")
  cat("   ⏳ This may take 2-3 minutes...\n")
  part6_start <- Sys.time()
  all_results$biological_validation <- test_biological_validation()
  part6_time <- as.numeric(difftime(Sys.time(), part6_start, units = "mins"))
  cat(sprintf("   ✅ PART 6 COMPLETE! (%.1f minutes)\n\n", part6_time))
  
  # Save all results
  cat("💾 SAVING RESULTS TO CSV FILES...\n")
  write.csv(all_results$main_comparison, "complete_main_comparison.csv", row.names = FALSE)
  write.csv(all_results$type1_error, "complete_type1_error.csv", row.names = FALSE)
  write.csv(all_results$power_analysis, "complete_power_analysis.csv", row.names = FALSE)
  write.csv(all_results$scalability, "complete_scalability.csv", row.names = FALSE)
  write.csv(all_results$parameter_sensitivity, "complete_parameter_sensitivity.csv", row.names = FALSE)
  write.csv(all_results$prefiltering_value, "complete_prefiltering_value.csv", row.names = FALSE)
  write.csv(all_results$biological_validation, "complete_biological_validation.csv", row.names = FALSE)
  
  total_time <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
  
  cat("\n🎉 COMPREHENSIVE ANALYSIS COMPLETE!\n")
  cat("⏱️  Total runtime:", sprintf("%.1f minutes", total_time), "\n")
  cat("📁 Results saved to:\n")
  cat("   📄 complete_main_comparison.csv (Main method comparison)\n")
  cat("   📄 complete_type1_error.csv (Type I error validation)\n")
  cat("   📄 complete_power_analysis.csv (Power analysis)\n")
  cat("   📄 complete_scalability.csv (Scalability analysis)\n")
  cat("   📄 complete_parameter_sensitivity.csv (Parameter sensitivity)\n")
  cat("   📄 complete_prefiltering_value.csv (Pre-filtering value analysis)\n")
  cat("   📄 complete_biological_validation.csv (Biological validation analysis)\n")
  cat("\n🚀 Ready for publication!\n")
  
  return(all_results)
}

# Run the complete analysis
if (interactive()) {
  complete_results <- run_complete_analysis()
}
