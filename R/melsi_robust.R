# -----------------------------------------------------------------------------
# MeLSI: ROBUST VERSION
# Simple, stable implementation that actually works
# -----------------------------------------------------------------------------

# Step 1: Install/Load necessary packages
if (!requireNamespace("vegan", quietly = TRUE)) install.packages("vegan")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")

library(vegan)
library(ggplot2)

# -----------------------------------------------------------------------------
# Function 1: Generate Synthetic Microbiome Data (UNCHANGED)
# -----------------------------------------------------------------------------
generate_test_data <- function(n_samples = 40, n_taxa = 100, n_signal_taxa = 8) {
    counts <- matrix(rnbinom(n_samples * n_taxa, mu = 50, size = 1), nrow = n_samples)
    colnames(counts) <- paste0("Taxon_", 1:n_taxa)
    rownames(counts) <- paste0("Sample_", 1:n_samples)
    
    metadata <- data.frame(
        SampleID = rownames(counts),
        Group = factor(rep(c("Control", "Treatment"), each = n_samples / 2))
    )
    
    # Inject signal
    treatment_indices <- which(metadata$Group == "Treatment")
    signal_taxa_indices <- 1:n_signal_taxa
    
    for (i in treatment_indices) {
        half <- floor(n_signal_taxa / 2)
        counts[i, 1:half] <- counts[i, 1:half] + rnbinom(half, mu = 25, size = 5)
        counts[i, (half + 1):n_signal_taxa] <- pmax(0, counts[i, (half + 1):n_signal_taxa] - rnbinom(n_signal_taxa - half, mu = 25, size = 5))
    }
    
    return(list(counts = counts, metadata = metadata))
}

# -----------------------------------------------------------------------------
# Function 2: Calculate PERMANOVA Pseudo-F Statistic (UNCHANGED)
# -----------------------------------------------------------------------------
calculate_permanova_F <- function(dist_matrix, labels) {
    permanova_res <- adonis2(dist_matrix ~ labels, permutations = 0)
    f_stat <- permanova_res$F[1]
    return(f_stat)
}

# -----------------------------------------------------------------------------
# Function 3: ROBUST Mahalanobis Distance
# -----------------------------------------------------------------------------
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

# -----------------------------------------------------------------------------
# Function 4: IMPROVED Weak Learner - Smarter and Faster
# -----------------------------------------------------------------------------
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

# -----------------------------------------------------------------------------
# Function 5: IMPROVED Ensemble Learner with Pre-filtering
# -----------------------------------------------------------------------------
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
        
        # Feature subsampling with stratification
        feature_indices <- sample(1:n_features, m, replace = FALSE)
        X_sub <- X[boot_indices, feature_indices, drop = FALSE]
        y_sub <- y[boot_indices]
        
        # Learn weak metric
        M_weak <- optimize_weak_learner_robust(X_sub, y_sub)
        
        # Calculate F-statistic for this weak learner
        dist_weak <- as.matrix(dist(X_sub %*% chol(M_weak)))
        f_stat_weak <- tryCatch({
            adonis2(dist_weak ~ y_sub, permutations = 0)$F[1]
        }, error = function(e) 0)
        
        # Only keep good weak learners
        if (f_stat_weak > 0) {
            M_full <- diag(ncol(X))  # Start with identity
            M_full[feature_indices, feature_indices] <- M_weak
            
            valid_count <- valid_count + 1
            learned_matrices[[valid_count]] <- M_full
            f_stats[valid_count] <- f_stat_weak
        }
        
        if (b %% 20 == 0) {
            cat(sprintf("Completed %d/%d weak learners (%.1f%%)\n", 
                       b, B, 100*b/B))
        }
    }
    
    if (valid_count == 0) return(diag(ncol(X)))
    
    # Weighted ensemble based on F-statistics
    weights <- f_stats[1:valid_count]
    weights <- weights / sum(weights)  # Normalize
    
    M_ensemble <- matrix(0, ncol(X), ncol(X))
    for (i in 1:valid_count) {
        M_ensemble <- M_ensemble + weights[i] * learned_matrices[[i]]
    }
    
    # Ensure positive definite
    eigen_M <- eigen(M_ensemble)
    eigen_M$values <- pmax(eigen_M$values, 0.01)
    M_ensemble <- eigen_M$vectors %*% diag(eigen_M$values) %*% t(eigen_M$vectors)
    
    cat(sprintf("Ensemble complete: %d valid learners, avg F-stat: %.4f\n", 
               valid_count, mean(f_stats[1:valid_count])))
    
    return(M_ensemble)
}

# -----------------------------------------------------------------------------
# Function 6: ROBUST Permutation Test
# -----------------------------------------------------------------------------
run_melsi_permutation_test <- function(X, y, n_perms = 19, B = 20, m_frac = 0.7, show_progress = TRUE) {
    if (show_progress) {
        cat("--- Starting MeLSI Analysis ---\n")
    }
    
    # 1. Learn metric on observed data
    if (show_progress) {
        cat("Learning metric on observed data...\n")
    }
    M_observed <- learn_melsi_metric_robust(X, y, B = B, m_frac = m_frac, pre_filter = TRUE)
    
    # Apply the same pre-filtering to X for distance calculation
    X_filtered <- X
    if (ncol(X) > 10) {
        # Replicate the pre-filtering logic
        feature_importance <- numeric(ncol(X))
        classes <- unique(y)
        if (length(classes) == 2) {
            class1_indices <- which(y == classes[1])
            class2_indices <- which(y == classes[2])
            
            for (i in 1:ncol(X)) {
                tryCatch({
                    test_result <- t.test(X[class1_indices, i], X[class2_indices, i])
                    feature_importance[i] <- abs(test_result$statistic)
                }, error = function(e) {
                    feature_importance[i] <- 0
                })
            }
            
            top_features <- order(feature_importance, decreasing = TRUE)
            n_keep <- max(10, min(ncol(X), floor(ncol(X) * 0.5)))
            keep_features <- top_features[1:n_keep]
            X_filtered <- X[, keep_features, drop = FALSE]
        }
    }
    
    dist_observed <- calculate_mahalanobis_dist_robust(X_filtered, M_observed)
    F_observed <- calculate_permanova_F(dist_observed, y)
    
    if (show_progress) {
        cat("Observed F-statistic:", F_observed, "\n")
    }
    
    # 2. Generate null distribution
    if (show_progress) {
        cat("Generating null distribution with", n_perms, "permutations...\n")
    }
    F_null <- numeric(n_perms)
    
    for (p in 1:n_perms) {
        y_permuted <- sample(y)
        M_permuted <- learn_melsi_metric_robust(X, y_permuted, B = B, m_frac = m_frac, pre_filter = TRUE)
        dist_permuted <- calculate_mahalanobis_dist_robust(X_filtered, M_permuted)
        F_null[p] <- calculate_permanova_F(dist_permuted, y_permuted)
        
        if (show_progress && p %% max(1, floor(n_perms/10)) == 0) {
            cat("Completed", p, "of", n_perms, "permutations\n")
        }
    }
    
    # 3. Calculate p-value
    p_value <- (sum(F_null >= F_observed) + 1) / (n_perms + 1)
    
    # 4. Diagnostics
    M_identity <- diag(ncol(X_filtered))
    metric_diff <- norm(M_observed - M_identity, 'F')
    eigen_vals <- eigen(M_observed)$values
    condition_number <- max(eigen_vals) / min(eigen_vals)
    
    if (show_progress) {
        cat("Analysis completed!\n")
        cat("Diagnostics:\n")
        cat("  Metric difference from identity:", round(metric_diff, 4), "\n")
        cat("  Condition number:", round(condition_number, 2), "\n")
        cat("  Eigenvalue range:", round(range(eigen_vals), 4), "\n")
    }
    
    return(list(
        F_observed = F_observed,
        p_value = p_value,
        M_ensemble = M_observed,
        null_distribution = F_null,
        diagnostics = list(
            metric_diff = metric_diff,
            condition_number = condition_number,
            eigenvalue_range = range(eigen_vals)
        )
    ))
}

# -----------------------------------------------------------------------------
# Function 7: Quick Test Function for Improved MeLSI
# -----------------------------------------------------------------------------
test_improved_melsi <- function() {
    cat("=== Testing Improved MeLSI ===\n")
    
    # Generate test data
    test_data <- generate_test_data(n_samples = 40, n_taxa = 50, n_signal_taxa = 8)
    X <- test_data$counts
    y <- test_data$metadata$Group
    
    # CLR transform
    X_clr <- X
    X_clr[X_clr == 0] <- 1e-10
    X_clr <- log(X_clr)
    X_clr <- X_clr - rowMeans(X_clr)
    
    # Test improved MeLSI
    cat("\n1. Testing Improved MeLSI...\n")
    start_time <- Sys.time()
    melsi_result <- run_melsi_permutation_test(X_clr, y, n_perms = 20, show_progress = TRUE)
    melsi_time <- Sys.time() - start_time
    
    # Test Euclidean baseline
    cat("\n2. Testing Euclidean baseline...\n")
    start_time <- Sys.time()
    dist_euclidean <- dist(X_clr)
    f_euclidean <- calculate_permanova_F(dist_euclidean, y)
    euclidean_time <- Sys.time() - start_time
    
    # Results
    cat("\n=== RESULTS ===\n")
    cat(sprintf("MeLSI F-statistic: %.4f (p = %.4f)\n", 
               melsi_result$F_observed, melsi_result$p_value))
    cat(sprintf("Euclidean F-statistic: %.4f\n", f_euclidean))
    cat(sprintf("MeLSI time: %.2f seconds\n", as.numeric(melsi_time)))
    cat(sprintf("Euclidean time: %.2f seconds\n", as.numeric(euclidean_time)))
    cat(sprintf("Speedup: %.1fx\n", as.numeric(melsi_time) / as.numeric(euclidean_time)))
    
    # Diagnostics
    cat(sprintf("\nMeLSI Diagnostics:\n"))
    cat(sprintf("  Metric difference from identity: %.4f\n", melsi_result$diagnostics$metric_diff))
    cat(sprintf("  Condition number: %.2f\n", melsi_result$diagnostics$condition_number))
    cat(sprintf("  Eigenvalue range: [%.4f, %.4f]\n", 
               melsi_result$diagnostics$eigenvalue_range[1], 
               melsi_result$diagnostics$eigenvalue_range[2]))
    
    return(list(
        melsi = melsi_result,
        euclidean = list(f_statistic = f_euclidean),
        times = list(melsi = melsi_time, euclidean = euclidean_time)
    ))
}

# Demo code removed - only run when explicitly called
# To run a demo, use: source('melsi_robust.R'); test_improved_melsi()
