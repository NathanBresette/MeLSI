# Improved MeLSI Algorithm - Fixed for Type I Error Control
# Key changes:
# 1. No pre-filtering bias (pre-filter on permuted data too)
# 2. Simpler ensemble to prevent overfitting
# 3. More conservative feature selection

source("melsi_robust.R")

# Improved MeLSI permutation test
run_melsi_improved <- function(X, y, n_perms = 99, B = 30, m_frac = 0.8, show_progress = TRUE) {
    if (show_progress) {
        cat("--- Starting IMPROVED MeLSI Analysis ---\n")
    }
    
    # 1. Learn metric on observed data (with conservative pre-filtering)
    if (show_progress) {
        cat("Learning metric on observed data...\n")
    }
    
    # Apply conservative pre-filtering (keep more features, less aggressive)
    X_filtered <- apply_conservative_prefiltering(X, y, filter_frac = 0.7)
    
    M_observed <- learn_melsi_metric_robust(X_filtered, y, B = B, m_frac = m_frac, pre_filter = FALSE)
    
    dist_observed <- calculate_mahalanobis_dist_robust(X_filtered, M_observed)
    F_observed <- calculate_permanova_F(dist_observed, y)
    
    if (show_progress) {
        cat("Observed F-statistic:", round(F_observed, 4), "\n")
    }
    
    # 2. Generate null distribution with CONSISTENT pre-filtering
    if (show_progress) {
        cat("Generating null distribution with", n_perms, "permutations...\n")
    }
    F_null <- numeric(n_perms)
    
    for (p in 1:n_perms) {
        # CRITICAL: Apply same pre-filtering to permuted data
        y_permuted <- sample(y)
        X_filtered_perm <- apply_conservative_prefiltering(X, y_permuted, filter_frac = 0.7)
        
        M_permuted <- learn_melsi_metric_robust(X_filtered_perm, y_permuted, B = B, m_frac = m_frac, pre_filter = FALSE)
        dist_permuted <- calculate_mahalanobis_dist_robust(X_filtered_perm, M_permuted)
        F_null[p] <- calculate_permanova_F(dist_permuted, y_permuted)
        
        if (show_progress && p %% max(1, floor(n_perms/10)) == 0) {
            cat("Completed", p, "of", n_perms, "permutations\n")
        }
    }
    
    # 3. Calculate p-value
    p_value <- (sum(F_null >= F_observed) + 1) / (n_perms + 1)
    
    if (show_progress) {
        cat("Analysis completed!\n")
        cat("F-statistic:", round(F_observed, 4), "\n")
        cat("P-value:", round(p_value, 4), "\n")
    }
    
    return(list(
        F_observed = F_observed,
        p_value = p_value,
        F_null = F_null,
        M_learned = M_observed,
        diagnostics = list(
            n_features_used = ncol(X_filtered),
            n_permutations = n_perms,
            ensemble_size = B
        )
    ))
}

# Conservative pre-filtering function
apply_conservative_prefiltering <- function(X, y, filter_frac = 0.7) {
    # Keep more features, use less aggressive filtering
    classes <- unique(y)
    if (length(classes) != 2 || ncol(X) <= 10) {
        return(X)
    }
    
    class1_indices <- which(y == classes[1])
    class2_indices <- which(y == classes[2])
    
    # Calculate feature importance with more conservative approach
    feature_importance <- numeric(ncol(X))
    for (i in 1:ncol(X)) {
        tryCatch({
            # Use variance instead of t-test to be less aggressive
            group1_var <- var(X[class1_indices, i])
            group2_var <- var(X[class2_indices, i])
            group1_mean <- mean(X[class1_indices, i])
            group2_mean <- mean(X[class2_indices, i])
            
            # Combine mean difference and variance for importance
            mean_diff <- abs(group1_mean - group2_mean)
            var_combined <- sqrt(group1_var + group2_var)
            feature_importance[i] <- mean_diff / (var_combined + 1e-10)
        }, error = function(e) {
            feature_importance[i] <- 0
        })
    }
    
    # Keep more features (70% instead of 50%)
    n_keep <- max(10, floor(ncol(X) * filter_frac))
    top_features <- order(feature_importance, decreasing = TRUE)[1:n_keep]
    
    return(X[, top_features, drop = FALSE])
}

# Test the improved algorithm
test_improved_melsi <- function() {
    cat("Testing IMPROVED MeLSI algorithm...\n\n")
    
    # Test 1: Null data (should NOT be significant)
    cat("=== TEST 1: NULL DATA (should NOT be significant) ===\n")
    set.seed(123)
    n_samples <- 100
    n_taxa <- 200
    
    # Create null data (no real group differences)
    counts_null <- matrix(rlnorm(n_samples * n_taxa, meanlog = 2, sdlog = 1), 
                         nrow = n_samples, ncol = n_taxa)
    metadata_null <- data.frame(
        Group = sample(c("Group1", "Group2"), n_samples, replace = TRUE)
    )
    
    # CLR transform
    counts_clr_null <- log(counts_null + 1)
    counts_clr_null <- counts_clr_null - rowMeans(counts_clr_null)
    
    # Test improved MeLSI
    results_null <- run_melsi_improved(
        counts_clr_null, metadata_null$Group,
        n_perms = 99, B = 30, m_frac = 0.8,
        show_progress = TRUE
    )
    
    cat("\nNULL DATA Results:\n")
    cat("F-statistic:", round(results_null$F_observed, 4), "\n")
    cat("P-value:", round(results_null$p_value, 4), "\n")
    cat("Significant:", ifelse(results_null$p_value < 0.05, "YES (PROBLEM!)", "NO (GOOD)"), "\n")
    
    # Test 2: Real data (should be significant)
    cat("\n=== TEST 2: REAL DATA (should be significant) ===\n")
    
    # Create data with real group differences
    counts_real <- counts_null
    group2_start <- n_samples/2 + 1
    counts_real[group2_start:n_samples, 1:30] <- 
        counts_real[group2_start:n_samples, 1:30] * 1.2  # 20% increase
    
    metadata_real <- data.frame(
        Group = c(rep("Group1", n_samples/2), rep("Group2", n_samples/2))
    )
    
    # CLR transform
    counts_clr_real <- log(counts_real + 1)
    counts_clr_real <- counts_clr_real - rowMeans(counts_clr_real)
    
    # Test improved MeLSI
    results_real <- run_melsi_improved(
        counts_clr_real, metadata_real$Group,
        n_perms = 99, B = 30, m_frac = 0.8,
        show_progress = TRUE
    )
    
    cat("\nREAL DATA Results:\n")
    cat("F-statistic:", round(results_real$F_observed, 4), "\n")
    cat("P-value:", round(results_real$p_value, 4), "\n")
    cat("Significant:", ifelse(results_real$p_value < 0.05, "YES (GOOD)", "NO (PROBLEM)"), "\n")
    
    # Summary
    cat("\n=== IMPROVED MeLSI SUMMARY ===\n")
    cat("Null data (no difference):", ifelse(results_null$p_value >= 0.05, "CORRECTLY NON-SIGNIFICANT", "FALSE POSITIVE"), "\n")
    cat("Real data (with difference):", ifelse(results_real$p_value < 0.05, "CORRECTLY SIGNIFICANT", "FALSE NEGATIVE"), "\n")
    
    success <- (results_null$p_value >= 0.05) && (results_real$p_value < 0.05)
    if (success) {
        cat("\n🎉 IMPROVED MeLSI WORKING CORRECTLY!\n")
    } else {
        cat("\n❌ IMPROVED MeLSI STILL HAS ISSUES\n")
    }
    
    return(list(
        null_results = results_null,
        real_results = results_real,
        success = success
    ))
}

# Run the test
if (interactive()) {
    improved_results <- test_improved_melsi()
}
