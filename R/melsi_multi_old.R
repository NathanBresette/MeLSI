# Multi-Group MeLSI Analysis Functions
# Provides both omnibus and pairwise analysis for 3+ groups

#' Run MeLSI Omnibus Analysis for Multiple Groups
#'
#' Performs omnibus MeLSI analysis to test for overall differences among 3+ groups.
#' Learns a global metric that maximizes separation across all group pairs.
#'
#' @param X A matrix of feature abundances with samples as rows and features as columns
#' @param y A vector of group labels for each sample
#' @param n_perms Number of permutations for p-value calculation (default: 75)
#' @param B Number of weak learners in the ensemble (default: 30)
#' @param m_frac Fraction of features to use in each weak learner (default: 0.8)
#' @param show_progress Whether to display progress information (default: TRUE)
#' @param plot_vip Whether to display a Variable Importance Plot (default: TRUE)
#'
#' @return A list containing:
#'   \item{F_observed}{The observed omnibus F-statistic}
#'   \item{p_value}{The permutation-based p-value}
#'   \item{feature_weights}{Global feature importance weights}
#'   \item{metric_matrix}{The learned global distance metric matrix}
#'   \item{group_info}{Information about groups analyzed}
#'
#' @export
melsi_omnibus <- function(X, y, n_perms = 75, B = 30, m_frac = 0.8, 
                         show_progress = TRUE, plot_vip = TRUE) {
    
    # Validate input
    groups <- unique(y)
    n_groups <- length(groups)
    
    if (n_groups < 3) {
        stop("melsi_omnibus requires 3 or more groups. Use melsi() for pairwise comparisons.")
    }
    
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
            cat("Permutation", p, "of", n_perms, "\n")
        }
    }
    
    # 5. Calculate p-value
    p_value <- (sum(F_null >= F_observed) + 1) / (n_perms + 1)
    
    # 6. Extract feature weights and display results
    feature_weights <- diag(M_observed)
    names(feature_weights) <- colnames(X_filtered)
    
    if (show_progress) {
        cat("Omnibus analysis completed!\n")
        cat("F-statistic:", round(F_observed, 4), "\n")
        cat("P-value:", round(p_value, 4), "\n")
        cat("\nGlobal feature importance: Access results$feature_weights to see which taxa\n")
        cat("contributed most to overall group separation.\n")
        
        # Show top 5 features
        if (length(feature_weights) >= 5) {
            top_5_idx <- order(feature_weights, decreasing = TRUE)[1:5]
            cat("\nTop 5 globally important features:\n")
            for (i in 1:5) {
                idx <- top_5_idx[i]
                feature_name <- if (!is.null(names(feature_weights)[idx])) {
                    names(feature_weights)[idx]
                } else {
                    paste0("Feature_", idx)
                }
                cat(sprintf("  %d. %s (weight: %.4f)\n", i, feature_name, feature_weights[idx]))
            }
        }
    }
    
    # 7. Create VIP plot if requested
    if (plot_vip && length(feature_weights) > 0) {
        tryCatch({
            plot_feature_importance(feature_weights, 
                                  main_title = paste0("Global Feature Importance (Omnibus)"))
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

#' Run MeLSI Pairwise Analysis for Multiple Groups
#'
#' Performs pairwise MeLSI analysis for all combinations of groups.
#' Provides comparison-specific feature importance for each pair.
#'
#' @param X A matrix of feature abundances with samples as rows and features as columns
#' @param y A vector of group labels for each sample
#' @param groups Optional vector of specific groups to compare (default: all groups)
#' @param n_perms Number of permutations for p-value calculation (default: 75)
#' @param B Number of weak learners in the ensemble (default: 30)
#' @param m_frac Fraction of features to use in each weak learner (default: 0.8)
#' @param show_progress Whether to display progress information (default: TRUE)
#' @param plot_vip Whether to display VIP plots for each comparison (default: FALSE)
#' @param correction_method Multiple testing correction method (default: "BH" for FDR)
#'
#' @return A list containing:
#'   \item{pairwise_results}{List of results for each pairwise comparison}
#'   \item{summary_table}{Summary table of all comparisons with corrected p-values}
#'   \item{significant_pairs}{List of significantly different pairs}
#'
#' @export
melsi_pairwise <- function(X, y, groups = NULL, n_perms = 75, B = 30, m_frac = 0.8,
                          show_progress = TRUE, plot_vip = FALSE, 
                          correction_method = "BH") {
    
    # Get all groups if not specified
    if (is.null(groups)) {
        groups <- unique(y)
    }
    
    n_groups <- length(groups)
    
    if (n_groups < 2) {
        stop("At least 2 groups required for pairwise analysis.")
    }
    
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
            melsi(X_pair, y_pair, n_perms = n_perms, B = B, m_frac = m_frac,
                  show_progress = FALSE, plot_vip = plot_vip)
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

#' Run Complete Multi-Group MeLSI Analysis
#'
#' Wrapper function that runs both omnibus and pairwise analyses for multiple groups.
#'
#' @param X A matrix of feature abundances with samples as rows and features as columns
#' @param y A vector of group labels for each sample
#' @param include_omnibus Whether to include omnibus analysis (default: TRUE)
#' @param include_pairwise Whether to include pairwise analysis (default: TRUE)
#' @param ... Additional arguments passed to omnibus and pairwise functions
#'
#' @return A list containing omnibus and pairwise results
#'
#' @export
melsi_multi <- function(X, y, include_omnibus = TRUE, include_pairwise = TRUE, ...) {
    
    groups <- unique(y)
    n_groups <- length(groups)
    
    if (n_groups < 2) {
        stop("At least 2 groups required. Use melsi() for single comparisons.")
    }
    
    cat("=== MeLSI Multi-Group Analysis ===\n")
    cat("Groups:", paste(groups, collapse = ", "), "\n")
    cat("Sample sizes:", paste(names(table(y)), ":", table(y), collapse = ", "), "\n\n")
    
    results <- list()
    
    # Run omnibus analysis if requested
    if (include_omnibus && n_groups >= 3) {
        cat("Running omnibus analysis...\n")
        # Filter out correction_method from omnibus parameters
        omnibus_params <- list(...)
        omnibus_params$correction_method <- NULL
        results$omnibus <- do.call(melsi_omnibus, c(list(X = X, y = y), omnibus_params))
    } else if (include_omnibus && n_groups == 2) {
        cat("Only 2 groups detected. Skipping omnibus analysis.\n")
    }
    
    # Run pairwise analysis if requested
    if (include_pairwise) {
        cat("\nRunning pairwise analysis...\n")
        results$pairwise <- do.call(melsi_pairwise, c(list(X = X, y = y), list(...)))
    }
    
    return(results)
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
    
    return(X[, top_features, drop = FALSE])
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
