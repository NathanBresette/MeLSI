# Consolidated MeLSI Analysis Function
# Handles both pairwise (2 groups) and multi-group (3+ groups) analysis

#' Run MeLSI Analysis
#'
#' Performs MeLSI (Metric Learning for Statistical Inference) analysis for microbiome data.
#' Automatically handles both pairwise comparisons (2 groups) and multi-group analysis (3+ groups).
#'
#' @param X A matrix of feature abundances with samples as rows and features as columns
#' @param y A vector of group labels for each sample
#' @param analysis_type Type of analysis to perform:
#'   - "auto" (default): Automatically choose based on number of groups
#'   - "pairwise": For 2 groups or all pairwise comparisons for 3+ groups
#'   - "omnibus": Global analysis for 3+ groups (requires at least 3 groups)
#'   - "both": Both omnibus and pairwise for 3+ groups
#' @param n_perms Number of permutations for p-value calculation (default: 75)
#' @param B Number of weak learners in the ensemble (default: 30)
#' @param m_frac Fraction of features to use in each weak learner (default: 0.8)
#' @param show_progress Whether to display progress information (default: TRUE)
#' @param plot_vip Whether to display Variable Importance Plot (default: TRUE)
#' @param correction_method Multiple testing correction method for pairwise comparisons (default: "BH")
#'
#' @return For 2 groups or pairwise analysis: List with F-statistic, p-value, feature weights, etc.
#'         For 3+ groups: List containing omnibus results, pairwise results, or both.
#'
#' @export
melsi <- function(X, y, analysis_type = "auto", n_perms = 75, B = 30, m_frac = 0.8, 
                 show_progress = TRUE, plot_vip = TRUE, correction_method = "BH") {
    
    # Validate input and ensure proper column names
    if (is.null(colnames(X)) || all(colnames(X) == "")) {
        colnames(X) <- paste0("Feature_", 1:ncol(X))
        if (show_progress) {
            cat("Warning: Input data has no column names. Using generic feature names.\n")
        }
    }
    
    groups <- unique(y)
    n_groups <- length(groups)
    
    if (n_groups < 2) {
        stop("MeLSI requires at least 2 groups.")
    }
    
    # Determine analysis type
    if (analysis_type == "auto") {
        if (n_groups == 2) {
            analysis_type <- "pairwise"
        } else {
            analysis_type <- "both"
        }
    }
    
    # Validate analysis type for number of groups
    if (n_groups == 2 && analysis_type %in% c("omnibus", "both")) {
        if (show_progress) {
            cat("Only 2 groups detected. Running pairwise analysis.\n")
        }
        analysis_type <- "pairwise"
    }
    
    if (n_groups >= 3 && analysis_type == "omnibus") {
        # Omnibus only for 3+ groups
        return(run_omnibus_analysis(X, y, n_perms, B, m_frac, show_progress, plot_vip))
    }
    
    if (analysis_type == "pairwise") {
        # Pairwise analysis
        if (n_groups == 2) {
            # Standard pairwise
            return(run_pairwise_analysis(X, y, n_perms, B, m_frac, show_progress, plot_vip))
        } else {
            # All pairwise comparisons
            return(run_all_pairwise_analysis(X, y, n_perms, B, m_frac, show_progress, 
                                           plot_vip, correction_method))
        }
    }
    
    if (analysis_type == "both") {
        # Both omnibus and pairwise
        if (show_progress) {
            cat("=== MeLSI Multi-Group Analysis ===\n")
            cat("Groups:", paste(groups, collapse = ", "), "\n")
            cat("Sample sizes:", paste(names(table(y)), ":", table(y), collapse = ", "), "\n\n")
        }
        
        results <- list()
        
        # Run omnibus analysis
        if (show_progress) cat("Running omnibus analysis...\n")
        results$omnibus <- run_omnibus_analysis(X, y, n_perms, B, m_frac, show_progress, plot_vip)
        
        # Run pairwise analysis
        if (show_progress) cat("\nRunning pairwise analysis...\n")
        results$pairwise <- run_all_pairwise_analysis(X, y, n_perms, B, m_frac, show_progress, 
                                                     plot_vip, correction_method)
        
        return(results)
    }
    
    stop("Invalid analysis_type. Use 'auto', 'pairwise', 'omnibus', or 'both'.")
}

# Helper function: Run standard pairwise analysis (2 groups)
run_pairwise_analysis <- function(X, y, n_perms, B, m_frac, show_progress, plot_vip) {
    
    if (show_progress) {
        cat("--- Starting MeLSI Analysis ---\n")
    }
    
    # 1. Learn metric on observed data (with conservative pre-filtering)
    if (show_progress) {
        cat("Learning metric on observed data...\n")
    }
    
    # Apply conservative pre-filtering
    X_filtered <- apply_conservative_prefiltering(X, y, filter_frac = 0.7)
    
    M_observed <- learn_melsi_metric_robust(X_filtered, y, B = B, m_frac = m_frac, pre_filter = FALSE)
    
    dist_observed <- calculate_mahalanobis_dist_robust(X_filtered, M_observed)
    F_observed <- calculate_permanova_F(dist_observed, y)
    
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
        
        if (show_progress) {
            cat("  [Permutation", p, "of", n_perms, "]\n")
            flush.console()
        }
    }
    
    # 3. Calculate p-value
    p_value <- (sum(F_null >= F_observed) + 1) / (n_perms + 1)
    
    # 4. Extract feature weights and calculate directionality
    feature_weights <- diag(M_observed)
    names(feature_weights) <- colnames(X_filtered)
    
    # Calculate directionality (which group has higher abundance)
    groups <- unique(y)
    directionality_info <- NULL
    mean_abundances <- NULL
    log2_fold_change <- NULL
    
    # Always calculate directionality for 2 groups
    if (length(groups) == 2) {
        group1_idx <- which(y == groups[1])
        group2_idx <- which(y == groups[2])
        
        # Ensure we have valid indices
        if (length(group1_idx) > 0 && length(group2_idx) > 0) {
            # Calculate mean abundances for each group
            mean_group1 <- colMeans(X_filtered[group1_idx, , drop = FALSE])
            mean_group2 <- colMeans(X_filtered[group2_idx, , drop = FALSE])
            
            # Determine directionality - ensure it's always a character vector with names
            directionality_info <- ifelse(mean_group1 > mean_group2, 
                                          paste0("Higher in ", as.character(groups[1])), 
                                          paste0("Higher in ", as.character(groups[2])))
            names(directionality_info) <- colnames(X_filtered)
            
            # Verify directionality was created correctly
            if (is.null(directionality_info) || length(directionality_info) == 0) {
                warning("Failed to calculate directionality. Creating default values.")
                directionality_info <- rep("Unknown", ncol(X_filtered))
                names(directionality_info) <- colnames(X_filtered)
            }
            
            # Calculate fold change and log2 fold change
            fold_change <- mean_group1 / (mean_group2 + 1e-10)
            log2_fold_change <- log2(fold_change)
            names(log2_fold_change) <- colnames(X_filtered)
            
            # Store mean abundances
            mean_abundances <- list(
                group1 = mean_group1,
                group2 = mean_group2,
                group1_name = as.character(groups[1]),
                group2_name = as.character(groups[2])
            )
        } else {
            warning("Invalid group indices. Cannot calculate directionality.")
        }
    }
    
    if (show_progress) {
        cat("Analysis completed!\n")
        cat("F-statistic:", round(F_observed, 4), "\n")
        cat("P-value:", round(p_value, 4), "\n")
        cat("\nFeature importance: Access results$feature_weights to see which taxa\n")
        cat("contributed most to the learned distance metric.\n")
        if (!is.null(directionality_info)) {
            cat("Directionality: Access results$directionality to see which group has higher abundance.\n")
        }
        
        # Show top 5 features with directionality if available
        if (length(feature_weights) >= 5) {
            top_5_idx <- order(feature_weights, decreasing = TRUE)[1:5]
            cat("\nTop 5 most important features:\n")
            for (i in 1:5) {
                idx <- top_5_idx[i]
                feature_name <- names(feature_weights)[idx]
                if (is.null(feature_name) || feature_name == "" || is.na(feature_name)) {
                    feature_name <- paste0("Feature_", idx)
                }
                direction_text <- ""
                if (!is.null(directionality_info) && !is.null(directionality_info[idx])) {
                    direction_text <- paste0(" [", directionality_info[idx], "]")
                }
                cat(sprintf("  %d. %s (weight: %.4f)%s\n", i, feature_name, feature_weights[idx], direction_text))
            }
        }
    }
    
    # 5. Create VIP plot if requested
    if (plot_vip && length(feature_weights) > 0) {
        tryCatch({
            plot_feature_importance(feature_weights, directionality = directionality_info)
        }, error = function(e) {
            if (show_progress) {
                cat("\nNote: Could not generate VIP plot. Error:", e$message, "\n")
            }
        })
    }
    
    # Return results - directionality should always be included for 2-group analysis
    # (will be NULL for multi-group, but should be a named vector for 2 groups)
    return(list(
        F_observed = F_observed,
        p_value = p_value,
        F_null = F_null,
        feature_weights = feature_weights,
        directionality = directionality_info,  # Named vector: "Higher in [group]" for each feature
        mean_abundances = mean_abundances,
        log2_fold_change = log2_fold_change,
        metric_matrix = M_observed,
        diagnostics = list(
            n_features_used = ncol(X_filtered),
            n_permutations = n_perms,
            ensemble_size = B
        )
    ))
}

# Helper function: Run omnibus analysis (3+ groups)
run_omnibus_analysis <- function(X, y, n_perms, B, m_frac, show_progress, plot_vip) {
    
    groups <- unique(y)
    n_groups <- length(groups)
    
    if (show_progress) {
        cat("--- Starting MeLSI Omnibus Analysis ---\n")
        cat("Groups:", paste(groups, collapse = ", "), "\n")
        cat("Sample sizes:", table(y), "\n")
    }
    
    # 1. Apply conservative pre-filtering
    if (show_progress) {
        cat("Applying conservative pre-filtering...\n")
    }
    
    X_filtered <- apply_conservative_prefiltering_multi(X, y, filter_frac = 0.7)
    
    # 2. Learn global metric optimized for all group pairs
    if (show_progress) {
        cat("Learning global metric for all group pairs...\n")
    }
    
    M_observed <- learn_melsi_metric_omnibus(X_filtered, y, B = B, m_frac = m_frac)
    
    # 3. Calculate omnibus F-statistic
    dist_observed <- calculate_mahalanobis_dist_robust(X_filtered, M_observed)
    F_observed <- calculate_permanova_F(dist_observed, y)
    
    # 4. Generate null distribution
    if (show_progress) {
        cat("Generating null distribution with", n_perms, "permutations...\n")
    }
    
    F_null <- numeric(n_perms)
    
    for (p in 1:n_perms) {
        y_permuted <- sample(y)
        X_filtered_perm <- apply_conservative_prefiltering_multi(X, y_permuted, filter_frac = 0.7)
        
        M_permuted <- learn_melsi_metric_omnibus(X_filtered_perm, y_permuted, B = B, m_frac = m_frac)
        dist_permuted <- calculate_mahalanobis_dist_robust(X_filtered_perm, M_permuted)
        F_null[p] <- calculate_permanova_F(dist_permuted, y_permuted)
        
        if (show_progress) {
            cat("  [Permutation", p, "of", n_perms, "]\n")
            flush.console()
        }
    }
    
    # 5. Calculate p-value
    p_value <- (sum(F_null >= F_observed) + 1) / (n_perms + 1)
    
    # 6. Extract feature weights and calculate directionality (highest mean group)
    feature_weights <- diag(M_observed)
    names(feature_weights) <- colnames(X_filtered)
    
    # Calculate which group has highest mean abundance for each feature
    directionality_info <- NULL
    mean_abundances <- NULL
    
    if (n_groups >= 2) {
        # Calculate mean abundances for each group
        mean_by_group <- list()
        for (g in groups) {
            group_idx <- which(y == g)
            mean_by_group[[as.character(g)]] <- colMeans(X_filtered[group_idx, , drop = FALSE])
        }
        
        # Determine which group has highest mean for each feature
        mean_matrix <- do.call(rbind, mean_by_group)
        max_group_idx <- apply(mean_matrix, 2, which.max)
        directionality_info <- as.character(groups[max_group_idx])
        names(directionality_info) <- colnames(X_filtered)
        
        # Store mean abundances
        mean_abundances <- mean_by_group
        names(mean_abundances) <- as.character(groups)
    }
    
    if (show_progress) {
        cat("Omnibus analysis completed!\n")
        cat("F-statistic:", round(F_observed, 4), "\n")
        cat("P-value:", round(p_value, 4), "\n")
        cat("\nGlobal feature importance: Access results$feature_weights to see which taxa\n")
        cat("contributed most to overall group separation.\n")
        if (!is.null(directionality_info)) {
            cat("Directionality: Access results$directionality to see which group has highest mean abundance.\n")
        }
        
        # Show top 5 features with directionality if available
        if (length(feature_weights) >= 5) {
            top_5_idx <- order(feature_weights, decreasing = TRUE)[1:5]
            cat("\nTop 5 globally important features:\n")
            for (i in 1:5) {
                idx <- top_5_idx[i]
                feature_name <- names(feature_weights)[idx]
                if (is.null(feature_name) || feature_name == "" || is.na(feature_name)) {
                    feature_name <- paste0("Feature_", idx)
                }
                direction_text <- ""
                if (!is.null(directionality_info) && !is.null(directionality_info[idx])) {
                    direction_text <- paste0(" [Highest in ", directionality_info[idx], "]")
                }
                cat(sprintf("  %d. %s (weight: %.4f)%s\n", i, feature_name, feature_weights[idx], direction_text))
            }
        }
    }
    
    # 7. Create VIP plot if requested
    if (plot_vip && length(feature_weights) > 0) {
        tryCatch({
            plot_feature_importance(feature_weights, 
                                  main_title = "Global Feature Importance (Omnibus)",
                                  directionality = directionality_info)
        }, error = function(e) {
            if (show_progress) {
                cat("\nNote: Could not generate VIP plot. Error:", e$message, "\n")
            }
        })
    }
    
    return(list(
        F_observed = F_observed,
        p_value = p_value,
        F_null = F_null,
        feature_weights = feature_weights,
        directionality = directionality_info,
        mean_abundances = mean_abundances,
        metric_matrix = M_observed,
        group_info = list(
            groups = groups,
            n_groups = n_groups,
            sample_sizes = table(y)
        ),
        diagnostics = list(
            n_features_used = ncol(X_filtered),
            n_permutations = n_perms,
            ensemble_size = B
        )
    ))
}

# Helper function: Run all pairwise comparisons (3+ groups)
run_all_pairwise_analysis <- function(X, y, n_perms, B, m_frac, show_progress, plot_vip, correction_method) {
    
    groups <- unique(y)
    n_groups <- length(groups)
    
    if (show_progress) {
        cat("--- Starting MeLSI Pairwise Analysis ---\n")
        cat("Groups:", paste(groups, collapse = ", "), "\n")
        cat("Total comparisons:", choose(n_groups, 2), "\n")
    }
    
    # Generate all pairwise combinations
    pairwise_combinations <- combn(groups, 2, simplify = FALSE)
    n_comparisons <- length(pairwise_combinations)
    
    # Store results
    pairwise_results <- list()
    summary_data <- data.frame(
        Group1 = character(n_comparisons),
        Group2 = character(n_comparisons),
        F_statistic = numeric(n_comparisons),
        P_value = numeric(n_comparisons),
        stringsAsFactors = FALSE
    )
    
    # Run pairwise comparisons
    for (i in 1:n_comparisons) {
        pair <- pairwise_combinations[[i]]
        group1 <- pair[1]
        group2 <- pair[2]
        
        if (show_progress) {
            cat(sprintf("\nComparison %d/%d: %s vs %s\n", i, n_comparisons, group1, group2))
        }
        
        # Subset data for this pair
        pair_indices <- y %in% pair
        X_pair <- X[pair_indices, , drop = FALSE]
        y_pair <- y[pair_indices]
        
        # Run MeLSI for this pair
        pair_result <- tryCatch({
            run_pairwise_analysis(X_pair, y_pair, n_perms = n_perms, B = B, m_frac = m_frac,
                                 show_progress = FALSE, plot_vip = FALSE)
        }, error = function(e) {
            if (show_progress) {
                cat("Error in", group1, "vs", group2, ":", e$message, "\n")
            }
            return(NULL)
        })
        
        if (!is.null(pair_result)) {
            # Store results
            comparison_name <- paste(group1, group2, sep = "_vs_")
            pairwise_results[[comparison_name]] <- pair_result
            
            # Update summary table
            summary_data$Group1[i] <- group1
            summary_data$Group2[i] <- group2
            summary_data$F_statistic[i] <- pair_result$F_observed
            summary_data$P_value[i] <- pair_result$p_value
            
            if (show_progress) {
                cat(sprintf("  F-statistic: %.4f, P-value: %.4f\n", 
                           pair_result$F_observed, pair_result$p_value))
            }
        }
    }
    
    # Apply multiple testing correction
    summary_data$P_value_corrected <- p.adjust(summary_data$P_value, method = correction_method)
    summary_data$Significant <- summary_data$P_value_corrected < 0.05
    
    # Identify significant pairs
    significant_pairs <- summary_data[summary_data$Significant, ]
    
    if (show_progress) {
        cat("\n=== Pairwise Analysis Summary ===\n")
        cat("Multiple testing correction:", correction_method, "\n")
        cat("Significant pairs (corrected p < 0.05):", nrow(significant_pairs), "\n")
        
        if (nrow(significant_pairs) > 0) {
            cat("\nSignificant comparisons:\n")
            for (i in 1:nrow(significant_pairs)) {
                row <- significant_pairs[i, ]
                cat(sprintf("  %s vs %s: F=%.3f, p=%.4f (corrected: %.4f)\n",
                           row$Group1, row$Group2, row$F_statistic, 
                           row$P_value, row$P_value_corrected))
            }
        }
    }
    
    return(list(
        pairwise_results = pairwise_results,
        summary_table = summary_data,
        significant_pairs = significant_pairs,
        correction_method = correction_method
    ))
}

# All helper functions from original melsi_robust.R
# (Conservative pre-filtering, metric learning, distance calculation, etc.)

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
    
    # Ensure column names are preserved
    filtered_X <- X[, top_features, drop = FALSE]
    colnames(filtered_X) <- colnames(X)[top_features]
    
    return(filtered_X)
}

# Helper function: Conservative pre-filtering for multi-group data
apply_conservative_prefiltering_multi <- function(X, y, filter_frac = 0.7) {
    classes <- unique(y)
    if (length(classes) < 2 || ncol(X) <= 10) {
        return(X)
    }
    
    # Calculate feature importance using ANOVA F-statistic for multi-group data
    feature_importance <- numeric(ncol(X))
    
    for (i in 1:ncol(X)) {
        tryCatch({
            # Use ANOVA F-statistic for multi-group comparison
            aov_result <- aov(X[, i] ~ y)
            f_stat <- summary(aov_result)[[1]][["F value"]][1]
            feature_importance[i] <- ifelse(is.na(f_stat), 0, f_stat)
        }, error = function(e) {
            feature_importance[i] <- 0
        })
    }
    
    # Keep top features
    n_keep <- max(10, floor(ncol(X) * filter_frac))
    top_features <- order(feature_importance, decreasing = TRUE)[1:n_keep]
    
    # Ensure column names are preserved
    filtered_X <- X[, top_features, drop = FALSE]
    colnames(filtered_X) <- colnames(X)[top_features]
    
    return(filtered_X)
}

# Helper function: Calculate PERMANOVA F-statistic
calculate_permanova_F <- function(dist_matrix, labels) {
    permanova_res <- vegan::adonis2(dist_matrix ~ labels, permutations = 0)
    f_stat <- permanova_res$F[1]
    return(f_stat)
}

# Helper function: Robust Mahalanobis Distance
calculate_mahalanobis_dist_robust <- function(X, M) {
    n_samples <- nrow(X)
    
    # Ensure M is positive definite
    eigen_M <- eigen(M)
    eigen_M$values <- pmax(eigen_M$values, 1e-6)  # Ensure positive eigenvalues
    
    # Compute M^(-1/2) safely
    M_half_inv <- eigen_M$vectors %*% diag(1/sqrt(eigen_M$values)) %*% t(eigen_M$vectors)
    
    # Transform data
    Y <- X %*% M_half_inv
    
    # Compute Euclidean distances
    dist_matrix <- as.matrix(dist(Y, method = "euclidean"))
    
    return(as.dist(dist_matrix))
}

# Helper function: Optimize weak learner
optimize_weak_learner_robust <- function(X, y, n_iterations = 50, learning_rate = 0.1) {
    n_samples <- nrow(X)
    n_features <- ncol(X)
    
    # Start with identity matrix
    M <- diag(n_features)
    
    # Get class information
    classes <- unique(y)
    if (length(classes) < 2) return(M)
    
    class1_indices <- which(y == classes[1])
    class2_indices <- which(y == classes[2])
    
    if (length(class1_indices) < 2 || length(class2_indices) < 2) return(M)
    
    # Track convergence for early stopping
    prev_f_stat <- -Inf
    stagnation_count <- 0
    max_stagnation <- 20
    
    # Simplified but improved gradient descent
    for (iter in 1:n_iterations) {
        # Sample one pair from each class (simpler but still effective)
        i1 <- sample(class1_indices, 1)
        j1 <- sample(setdiff(class1_indices, i1), 1)
        i2 <- sample(class2_indices, 1)
        j2 <- sample(setdiff(class2_indices, i2), 1)
        
        # Compute differences
        diff1 <- X[i1, ] - X[j1, ]  # Within class 1
        diff2 <- X[i2, ] - X[j2, ]  # Within class 2
        diff3 <- X[i1, ] - X[i2, ]  # Between classes
        
        # Vectorized gradient calculation (the key improvement!)
        grad_between <- diff3^2
        grad_within <- -(diff1^2 + diff2^2) / 2
        total_gradient <- grad_between + grad_within
        
        # Adaptive learning rate
        current_learning_rate <- learning_rate * (1 / (1 + iter * 0.1))
        
        # Vectorized update (much faster than the old for loop!)
        diag(M) <- diag(M) + current_learning_rate * total_gradient
        diag(M) <- pmax(diag(M), 0.01)  # Keep positive
        
        # Early stopping if no improvement
        if (iter %% 20 == 0) {
            dist_matrix <- as.matrix(dist(X %*% chol(M)))
            current_f_stat <- tryCatch({
                vegan::adonis2(dist_matrix ~ y, permutations = 0)$F[1]
            }, error = function(e) 0)
            
            if (current_f_stat <= prev_f_stat) {
                stagnation_count <- stagnation_count + 1
                if (stagnation_count >= 5) break
            } else {
                stagnation_count <- 0
            }
            prev_f_stat <- current_f_stat
        }
    }
    
    return(M)
}

# Helper function: Learn MeLSI metric
learn_melsi_metric_robust <- function(X, y, B = 20, m_frac = 0.7, 
                                     pre_filter = TRUE, 
                                     filter_threshold = 0.1) {
    n_samples <- nrow(X)
    n_features <- ncol(X)
    m <- max(2, floor(n_features * m_frac))
    
    # Pre-filtering: Remove features with low variance or no signal
    if (pre_filter && n_features > 10) {
        # Calculate feature importance using simple t-test
        feature_importance <- numeric(n_features)
        classes <- unique(y)
        if (length(classes) == 2) {
            class1_indices <- which(y == classes[1])
            class2_indices <- which(y == classes[2])
            
            for (i in 1:n_features) {
                # Simple t-test for feature importance
                tryCatch({
                    test_result <- t.test(X[class1_indices, i], X[class2_indices, i])
                    feature_importance[i] <- abs(test_result$statistic)
                }, error = function(e) {
                    feature_importance[i] <- 0
                })
            }
            
            # Keep top features
            top_features <- order(feature_importance, decreasing = TRUE)
            n_keep <- max(10, min(n_features, floor(n_features * 0.5)))
            keep_features <- top_features[1:n_keep]
            
            X <- X[, keep_features, drop = FALSE]
            n_features <- ncol(X)
            m <- max(2, floor(n_features * m_frac))
        }
    }
    
    learned_matrices <- vector("list", B)
    valid_count <- 0
    f_stats <- numeric(B)
    
    for (b in 1:B) {
        # Bootstrap sampling
        boot_indices <- sample(1:n_samples, n_samples, replace = TRUE)
        if (length(unique(y[boot_indices])) < 2) next
        
        X_boot <- X[boot_indices, , drop = FALSE]
        y_boot <- y[boot_indices]
        
        # Feature subsampling
        feature_indices <- sample(1:n_features, m, replace = FALSE)
        X_subset <- X_boot[, feature_indices, drop = FALSE]
        
        # Learn weak learner
        M_weak <- optimize_weak_learner_robust(X_subset, y_boot)
        
        # Expand back to full feature space
        M_full <- diag(n_features)
        M_full[feature_indices, feature_indices] <- M_weak
        
        # Validate weak learner
        tryCatch({
            dist_test <- calculate_mahalanobis_dist_robust(X_subset, M_weak)
            f_stat <- calculate_permanova_F(dist_test, y_boot)
            
            if (is.finite(f_stat) && f_stat > 0) {
                valid_count <- valid_count + 1
                learned_matrices[[valid_count]] <- M_full
                f_stats[valid_count] <- f_stat
            }
        }, error = function(e) {
            # Skip this weak learner if it fails
        })
        
        if (valid_count >= B) break
    }
    
    if (valid_count == 0) {
        warning("No valid weak learners found. Returning identity matrix.")
        return(diag(n_features))
    }
    
    # Ensemble averaging with performance weighting
    weights <- f_stats[1:valid_count]
    weights <- weights / sum(weights)  # Normalize weights
    
    M_ensemble <- matrix(0, n_features, n_features)
    for (i in 1:valid_count) {
        M_ensemble <- M_ensemble + weights[i] * learned_matrices[[i]]
    }
    
    # Ensure positive definiteness
    eigen_result <- eigen(M_ensemble)
    eigen_result$values <- pmax(eigen_result$values, 1e-6)
    M_ensemble <- eigen_result$vectors %*% diag(eigen_result$values) %*% t(eigen_result$vectors)
    
    return(M_ensemble)
}

# Helper function: Learn omnibus metric for multi-group data
learn_melsi_metric_omnibus <- function(X, y, B = 30, m_frac = 0.8) {
    n_samples <- nrow(X)
    n_features <- ncol(X)
    m <- max(2, floor(n_features * m_frac))
    
    learned_matrices <- vector("list", B)
    valid_count <- 0
    f_stats <- numeric(B)
    
    for (b in 1:B) {
        # Bootstrap sampling
        boot_indices <- sample(1:n_samples, n_samples, replace = TRUE)
        if (length(unique(y[boot_indices])) < 2) next
        
        X_boot <- X[boot_indices, , drop = FALSE]
        y_boot <- y[boot_indices]
        
        # Feature subsampling
        feature_indices <- sample(1:n_features, m, replace = FALSE)
        X_subset <- X_boot[, feature_indices, drop = FALSE]
        
        # Learn weak learner optimized for all group pairs
        M_weak <- optimize_weak_learner_omnibus(X_subset, y_boot)
        
        # Expand back to full feature space
        M_full <- diag(n_features)
        M_full[feature_indices, feature_indices] <- M_weak
        
        # Validate weak learner
        tryCatch({
            dist_test <- calculate_mahalanobis_dist_robust(X_subset, M_weak)
            f_stat <- calculate_permanova_F(dist_test, y_boot)
            
            if (is.finite(f_stat) && f_stat > 0) {
                valid_count <- valid_count + 1
                learned_matrices[[valid_count]] <- M_full
                f_stats[valid_count] <- f_stat
            }
        }, error = function(e) {
            # Skip this weak learner if it fails
        })
        
        if (valid_count >= B) break
    }
    
    if (valid_count == 0) {
        warning("No valid weak learners found. Returning identity matrix.")
        return(diag(n_features))
    }
    
    # Ensemble averaging with performance weighting
    weights <- f_stats[1:valid_count]
    weights <- weights / sum(weights)
    
    M_ensemble <- matrix(0, n_features, n_features)
    for (i in 1:valid_count) {
        M_ensemble <- M_ensemble + weights[i] * learned_matrices[[i]]
    }
    
    # Ensure positive definiteness
    eigen_result <- eigen(M_ensemble)
    eigen_result$values <- pmax(eigen_result$values, 1e-6)
    M_ensemble <- eigen_result$vectors %*% diag(eigen_result$values) %*% t(eigen_result$vectors)
    
    return(M_ensemble)
}

# Helper function: Optimize weak learner for omnibus analysis
optimize_weak_learner_omnibus <- function(X, y, n_iterations = 50, learning_rate = 0.1) {
    n_samples <- nrow(X)
    n_features <- ncol(X)
    
    # Start with identity matrix
    M <- diag(n_features)
    
    # Get all group pairs
    groups <- unique(y)
    n_groups <- length(groups)
    
    if (n_groups < 2) return(M)
    
    # Track convergence
    prev_f_stat <- -Inf
    stagnation_count <- 0
    
    for (iter in 1:n_iterations) {
        # Sample from all group pairs (balanced sampling)
        group_pairs <- combn(groups, 2, simplify = FALSE)
        
        # Randomly select a group pair to optimize for this iteration
        selected_pair <- sample(group_pairs, 1)[[1]]
        group1_indices <- which(y == selected_pair[1])
        group2_indices <- which(y == selected_pair[2])
        
        if (length(group1_indices) < 2 || length(group2_indices) < 2) next
        
        # Sample from selected pair (same as pairwise optimization)
        i1 <- sample(group1_indices, 1)
        j1 <- sample(setdiff(group1_indices, i1), 1)
        i2 <- sample(group2_indices, 1)
        j2 <- sample(setdiff(group2_indices, i2), 1)
        
        # Compute differences
        diff1 <- X[i1, ] - X[j1, ]  # Within group 1
        diff2 <- X[i2, ] - X[j2, ]  # Within group 2
        diff3 <- X[i1, ] - X[i2, ]  # Between groups
        
        # Gradient calculation
        grad_between <- diff3^2
        grad_within <- -(diff1^2 + diff2^2) / 2
        total_gradient <- grad_between + grad_within
        
        # Adaptive learning rate
        current_learning_rate <- learning_rate * (1 / (1 + iter * 0.1))
        
        # Update
        diag(M) <- diag(M) + current_learning_rate * total_gradient
        diag(M) <- pmax(diag(M), 0.01)
        
        # Early stopping check
        if (iter %% 20 == 0) {
            dist_matrix <- as.matrix(dist(X %*% chol(M)))
            current_f_stat <- tryCatch({
                vegan::adonis2(dist_matrix ~ y, permutations = 0)$F[1]
            }, error = function(e) 0)
            
            if (current_f_stat <= prev_f_stat) {
                stagnation_count <- stagnation_count + 1
                if (stagnation_count >= 5) break
            } else {
                stagnation_count <- 0
            }
            prev_f_stat <- current_f_stat
        }
    }
    
    return(M)
}

#' Plot Feature Importance from MeLSI Analysis
#'
#' Creates a barplot showing the top features ranked by their learned weights
#'
#' @param feature_weights Named vector of feature weights
#' @param top_n Number of top features to display (default: 8)
#' @param main_title Optional title for the plot
#' @param directionality Optional named vector indicating which group has higher abundance for each feature
#'
#' @export
plot_feature_importance <- function(feature_weights, top_n = 8, main_title = NULL, directionality = NULL) {
    # Validate input
    if (length(feature_weights) == 0) {
        stop("feature_weights is empty")
    }
    
    # Load required packages
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("ggplot2 package is required for plotting. Please install it with: install.packages('ggplot2')")
    }
    
    # Sort features by weight
    sorted_weights <- sort(feature_weights, decreasing = TRUE)
    
    # Take top N features
    n_display <- min(top_n, length(sorted_weights))
    top_weights <- sorted_weights[1:n_display]
    
    # Get feature names
    feature_names <- names(top_weights)
    if (is.null(feature_names)) {
        feature_names <- paste0("Feature_", 1:n_display)
    }
    
    # Truncate long names for better display (more aggressive truncation)
    feature_names <- ifelse(nchar(feature_names) > 20, 
                           paste0(substr(feature_names, 1, 17), "..."),
                           feature_names)
    
    # Extract directionality for top features if provided
    directionality_colors <- NULL
    directionality_labels <- NULL
    if (!is.null(directionality) && length(directionality) > 0) {
        # Match directionality to top features
        directionality_labels <- directionality[names(top_weights)]
        
        # Create color mapping based on directionality
        # Extract unique group names from directionality strings
        unique_groups <- unique(directionality_labels)
        if (length(unique_groups) == 2) {
            # Two groups - use two colors
            directionality_colors <- ifelse(directionality_labels == unique_groups[1], 
                                           "#E63946",  # Red for group 1
                                           "#457B9D")  # Blue for group 2
        } else {
            # More than 2 groups or missing - use default color
            directionality_colors <- rep("steelblue", length(directionality_labels))
        }
    }
    
    # Create data frame for ggplot
    plot_data <- data.frame(
        Feature = factor(feature_names, levels = rev(feature_names)),  # Reverse for top-to-bottom ordering
        Weight = top_weights
    )
    
    # Add directionality info if available
    if (!is.null(directionality_labels)) {
        plot_data$Directionality <- directionality_labels
        plot_data$Color <- directionality_colors
    }
    
    # Create title
    title_text <- if (!is.null(main_title)) main_title else paste0("Top ", n_display, " Features by Importance")
    
    # Create ggplot with or without directionality coloring
    if (!is.null(directionality_labels) && !is.null(directionality_colors)) {
        # Plot with directionality colors
        p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Weight, y = Feature, fill = Directionality)) +
            ggplot2::geom_col() +
            ggplot2::scale_fill_manual(values = setNames(directionality_colors, unique(directionality_labels)),
                                      name = "Higher in") +
            ggplot2::geom_text(ggplot2::aes(label = sprintf("%.3f", Weight)), 
                              hjust = -0.1, size = 3) +
            ggplot2::labs(
                title = title_text,
                x = "Feature Weight",
                y = ""
            ) +
            ggplot2::theme_bw() +
            ggplot2::theme(
                plot.title = ggplot2::element_text(size = 12, hjust = 0.5),
                axis.text.y = ggplot2::element_text(size = 10),
                axis.text.x = ggplot2::element_text(size = 9),
                axis.title.x = ggplot2::element_text(size = 11),
                panel.grid.minor = ggplot2::element_blank(),
                plot.margin = ggplot2::margin(10, 30, 10, 10),
                legend.position = "right"
            ) +
            ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.1)))
    } else {
        # Plot without directionality (original behavior)
        p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Weight, y = Feature)) +
            ggplot2::geom_col(fill = "steelblue") +
            ggplot2::geom_text(ggplot2::aes(label = sprintf("%.3f", Weight)), 
                              hjust = -0.1, size = 3) +
            ggplot2::labs(
                title = title_text,
                x = "Feature Weight",
                y = ""
            ) +
            ggplot2::theme_bw() +
            ggplot2::theme(
                plot.title = ggplot2::element_text(size = 12, hjust = 0.5),
                axis.text.y = ggplot2::element_text(size = 10),
                axis.text.x = ggplot2::element_text(size = 9),
                axis.title.x = ggplot2::element_text(size = 11),
                panel.grid.minor = ggplot2::element_blank(),
                plot.margin = ggplot2::margin(10, 30, 10, 10)
            ) +
            ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.1)))
    }
    
    # Print the plot
    print(p)
    
    # Return the plot object for further customization if needed
    return(invisible(p))
}
