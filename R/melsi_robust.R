# Improved MeLSI Algorithm - Fixed for Type I Error Control
# Key changes:
# 1. No pre-filtering bias (pre-filter on permuted data too)
# 2. Simpler ensemble to prevent overfitting
# 3. More conservative feature selection

# Load required packages
if (!requireNamespace("vegan", quietly = TRUE)) install.packages("vegan")
library(vegan)

#' Run MeLSI Analysis
#'
#' Performs MeLSI (Metric Learning for Statistical Inference) analysis with permutation testing
#' for microbiome data.
#'
#' @param X A matrix of feature abundances with samples as rows and features as columns
#' @param y A vector of group labels for each sample
#' @param n_perms Number of permutations for p-value calculation (default: 75)
#' @param B Number of weak learners in the ensemble (default: 30)
#' @param m_frac Fraction of features to use in each weak learner (default: 0.8)
#' @param show_progress Whether to display progress information (default: TRUE)
#'
#' @return A list containing:
#'   \item{F_observed}{The observed F-statistic}
#'   \item{p_value}{The permutation-based p-value}
#'   \item{feature_weights}{The learned feature weights}
#'   \item{metric_matrix}{The learned distance metric matrix}
#'
#' @export
melsi <- function(X, y, n_perms = 75, B = 30, m_frac = 0.8, show_progress = TRUE) {
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
        feature_weights = diag(M_observed),
        metric_matrix = M_observed,
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

# Helper function: Calculate PERMANOVA F-statistic
calculate_permanova_F <- function(dist_matrix, labels) {
    permanova_res <- adonis2(dist_matrix ~ labels, permutations = 0)
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
                adonis2(dist_matrix ~ y, permutations = 0)$F[1]
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
            
            cat(sprintf("Pre-filtering: Keeping %d/%d features (%.1f%%)\n", 
                       n_keep, n_features, 100*n_keep/n_features))
            
            X <- X[, keep_features, drop = FALSE]
            n_features <- ncol(X)
            m <- max(2, floor(n_features * m_frac))
        }
    }
    
    learned_matrices <- vector("list", B)
    valid_count <- 0
    f_stats <- numeric(B)
    
    cat(sprintf("Training %d weak learners with %d features each...\n", B, m))
    
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
        cat("Warning: No valid weak learners found. Returning identity matrix.\n")
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
    
    cat(sprintf("Successfully trained %d weak learners\n", valid_count))
    return(M_ensemble)
}