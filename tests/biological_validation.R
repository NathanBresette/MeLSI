# Biological Validation Analysis for MeLSI
# This script analyzes learned metric weights for biological relevance

# Load required packages
library(phyloseq)
library(vegan)
library(microbiome)

# Source MeLSI functions
source("melsi_robust.R")

# Function: Biological Validation Analysis
analyze_biological_validation <- function() {
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
    
    # Extract counts and metadata
    counts <- as.matrix(otu_table(sex_subset))
    counts <- t(counts)  # Transpose to samples x taxa
    
    # Ensure counts is numeric matrix and has proper dimensions
    counts <- matrix(as.numeric(counts), nrow = nrow(counts), ncol = ncol(counts))
    rownames(counts) <- sample_names(sex_subset)
    colnames(counts) <- taxa_names(sex_subset)
    
    metadata <- data.frame(
      SampleID = sample_names(sex_subset),
      Group = sample_data(sex_subset)$sex
    )
    
    # CLR transformation
    counts_clr <- counts
    counts_clr[counts_clr == 0] <- 1e-10
    counts_clr <- log(counts_clr)
    counts_clr <- counts_clr - rowMeans(counts_clr, na.rm = TRUE)
    
    # Ensure no NaN or Inf values
    counts_clr[is.na(counts_clr)] <- 0
    counts_clr[is.infinite(counts_clr)] <- 0
    
    # Run MeLSI
    cat("  Running MeLSI on Atlas1006 sex data...\n")
    cat("  Data dimensions:", nrow(counts_clr), "samples x", ncol(counts_clr), "taxa\n")
    cat("  Group levels:", paste(unique(metadata$Group), collapse = ", "), "\n")
    
    # Check for any issues with the data
    if (any(is.na(counts_clr)) || any(is.infinite(counts_clr))) {
      cat("  Warning: Found NA or Inf values in CLR data\n")
    }
    
    melsi_results <- run_melsi_permutation_test(
      counts_clr, metadata$Group, 
      n_perms = 19, B = 20, m_frac = 0.7, 
      show_progress = FALSE
    )
    
    # Extract learned metric weights
    M_learned <- melsi_results$M_ensemble
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
    
    # Check that we have the right number of groups
    unique_groups <- unique(metadata$Group)
    if (length(unique_groups) != 2) {
      cat("  Warning: Expected 2 groups, found", length(unique_groups), "\n")
      stop("Need exactly 2 groups for comparison")
    }
    
    group1_idx <- metadata$Group == unique_groups[1]
    group2_idx <- metadata$Group == unique_groups[2]
    
    # Check that we have samples in both groups
    if (sum(group1_idx) == 0 || sum(group2_idx) == 0) {
      cat("  Error: One or both groups have no samples\n")
      stop("Both groups must have samples")
    }
    
    cat("  Group 1 (", unique_groups[1], "):", sum(group1_idx), "samples\n")
    cat("  Group 2 (", unique_groups[2], "):", sum(group2_idx), "samples\n")
    
    # Simple t-test for each taxon
    p_values <- numeric(ncol(counts_clr))
    effect_sizes <- numeric(ncol(counts_clr))
    
    for (i in 1:ncol(counts_clr)) {
      tryCatch({
        test_result <- t.test(counts_clr[group1_idx, i], counts_clr[group2_idx, i])
        p_values[i] <- test_result$p.value
        effect_sizes[i] <- abs(test_result$statistic)
      }, error = function(e) {
        p_values[i] <<- 1
        effect_sizes[i] <<- 0
      })
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
    tryCatch({
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
      cat("  Error storing results:", e$message, "\n")
      cat("  Dimensions - samples:", nrow(counts_clr), "taxa:", ncol(counts_clr), "\n")
      cat("  Metric weights length:", length(metric_weights), "\n")
      stop("Failed to store results")
    })
    
    cat("✅ Atlas1006 analysis completed successfully!\n")
    
  }, error = function(e) {
    cat("❌ Could not analyze Atlas1006 sex data:", e$message, "\n")
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
    
    # Extract counts and metadata
    counts <- as.matrix(otu_table(combined_subset))
    counts <- t(counts)  # Transpose to samples x taxa
    
    # Ensure counts is numeric matrix and has proper dimensions
    counts <- matrix(as.numeric(counts), nrow = nrow(counts), ncol = ncol(counts))
    rownames(counts) <- sample_names(combined_subset)
    colnames(counts) <- taxa_names(combined_subset)
    
    metadata <- data.frame(
      SampleID = sample_names(combined_subset),
      Group = sample_data(combined_subset)$warmed
    )
    
    # CLR transformation
    counts_clr <- counts
    counts_clr[counts_clr == 0] <- 1e-10
    counts_clr <- log(counts_clr)
    counts_clr <- counts_clr - rowMeans(counts_clr, na.rm = TRUE)
    
    # Ensure no NaN or Inf values
    counts_clr[is.na(counts_clr)] <- 0
    counts_clr[is.infinite(counts_clr)] <- 0
    
    # Run MeLSI
    cat("  Running MeLSI on SoilRep warming data...\n")
    cat("  Data dimensions:", nrow(counts_clr), "samples x", ncol(counts_clr), "taxa\n")
    cat("  Group levels:", paste(unique(metadata$Group), collapse = ", "), "\n")
    
    # Check for any issues with the data
    if (any(is.na(counts_clr)) || any(is.infinite(counts_clr))) {
      cat("  Warning: Found NA or Inf values in CLR data\n")
    }
    
    melsi_results <- run_melsi_permutation_test(
      counts_clr, metadata$Group, 
      n_perms = 19, B = 20, m_frac = 0.7, 
      show_progress = FALSE
    )
    
    # Extract learned metric weights
    M_learned <- melsi_results$M_ensemble
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
    
    cat("✅ SoilRep analysis completed successfully!\n")
    
  }, error = function(e) {
    cat("❌ Could not analyze SoilRep warming data:", e$message, "\n")
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

# Run the analysis
cat("🚀 Starting Biological Validation Analysis...\n")
biological_results <- analyze_biological_validation()

# Save results
write.csv(biological_results, "biological_validation_results.csv", row.names = FALSE)
cat("\n💾 Results saved to biological_validation_results.csv\n")

cat("\n🎉 Biological Validation Analysis Complete!\n")
