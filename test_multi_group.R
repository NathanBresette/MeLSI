# Test script for multi-group MeLSI analysis
# Creates synthetic data with 3+ groups and tests omnibus vs pairwise approaches

# Load required packages
library(MeLSI)

# Set seed for reproducibility
set.seed(42)

# Function to create synthetic multi-group microbiome data
create_multi_group_data <- function(n_samples_per_group = 40, n_taxa = 200, 
                                   signal_strength = "medium", noise_level = 0.1) {
    
    n_groups <- 4  # Control, Mild, Moderate, Severe
    total_samples <- n_samples_per_group * n_groups
    group_names <- c("Control", "Mild", "Moderate", "Severe")
    
    cat("Creating synthetic data with", n_groups, "groups and", n_taxa, "taxa\n")
    cat("Groups:", paste(group_names, collapse = ", "), "\n")
    cat("Samples per group:", n_samples_per_group, "\n")
    
    # Create base microbiome data (log-normal distribution)
    X <- matrix(rlnorm(total_samples * n_taxa, meanlog = 2, sdlog = 1), 
                nrow = total_samples, ncol = n_taxa)
    
    # Add realistic noise
    X <- X + matrix(rnorm(total_samples * n_taxa, 0, noise_level), 
                    nrow = total_samples, ncol = n_taxa)
    X <- pmax(X, 0)  # Ensure non-negative
    
    # Create group labels
    y <- rep(group_names, each = n_samples_per_group)
    
    # Define signal strength
    if (signal_strength == "weak") {
        multiplier <- c(1.0, 1.2, 1.3, 1.4)  # Subtle differences
        taxa_affected <- 15  # Few taxa affected
    } else if (signal_strength == "medium") {
        multiplier <- c(1.0, 1.4, 1.7, 2.0)  # Moderate differences
        taxa_affected <- 25  # Moderate number of taxa
    } else if (signal_strength == "strong") {
        multiplier <- c(1.0, 1.8, 2.5, 3.0)  # Strong differences
        taxa_affected <- 35  # Many taxa affected
    }
    
    # Add progressive signal across groups
    # Different taxa respond to different group comparisons
    
    # Taxa 1-10: Progressive increase across all groups (global signal)
    for (i in 1:min(10, taxa_affected)) {
        for (g in 1:n_groups) {
            group_start <- (g-1) * n_samples_per_group + 1
            group_end <- g * n_samples_per_group
            X[group_start:group_end, i] <- X[group_start:group_end, i] * multiplier[g]
        }
    }
    
    # Taxa 11-20: Control vs Others (pairwise signal)
    for (i in 11:min(20, taxa_affected)) {
        # Control group stays baseline
        # Other groups increase
        for (g in 2:n_groups) {
            group_start <- (g-1) * n_samples_per_group + 1
            group_end <- g * n_samples_per_group
            X[group_start:group_end, i] <- X[group_start:group_end, i] * (1.5 + (g-2) * 0.3)
        }
    }
    
    # Taxa 21-30: Mild vs Severe (specific pairwise signal)
    if (taxa_affected >= 25) {
        for (i in 21:min(30, taxa_affected)) {
            # Mild group
            mild_start <- n_samples_per_group + 1
            mild_end <- 2 * n_samples_per_group
            X[mild_start:mild_end, i] <- X[mild_start:mild_end, i] * 1.3
            
            # Severe group  
            severe_start <- 3 * n_samples_per_group + 1
            severe_end <- 4 * n_samples_per_group
            X[severe_start:severe_end, i] <- X[severe_start:severe_end, i] * 2.2
        }
    }
    
    # CLR transformation (recommended for microbiome data)
    X_clr <- log(X + 1)
    X_clr <- X_clr - rowMeans(X_clr)
    
    # Add column names for better interpretability
    colnames(X_clr) <- paste0("Taxon_", sprintf("%03d", 1:n_taxa))
    
    cat("Signal added to", taxa_affected, "taxa\n")
    cat("Signal strength:", signal_strength, "\n")
    
    return(list(
        X = X_clr,
        y = y,
        group_names = group_names,
        taxa_affected = taxa_affected,
        signal_info = list(
            strength = signal_strength,
            multiplier = multiplier,
            global_taxa = 1:min(10, taxa_affected),
            control_vs_others_taxa = 11:min(20, taxa_affected),
            mild_vs_severe_taxa = 21:min(30, taxa_affected)
        )
    ))
}

# Function to run comprehensive multi-group test
run_multi_group_test <- function(signal_strength = "medium") {
    
    cat("=== MeLSI Multi-Group Test ===\n")
    cat("Signal strength:", signal_strength, "\n\n")
    
    # Create synthetic data
    data <- create_multi_group_data(signal_strength = signal_strength)
    X <- data$X
    y <- data$y
    
    cat("\n=== Running MeLSI Multi-Group Analysis ===\n")
    
    # Run complete multi-group analysis
    results <- melsi_multi(
        X = X, 
        y = y,
        include_omnibus = TRUE,
        include_pairwise = TRUE,
        n_perms = 50,  # Reduced for faster testing
        B = 20,         # Reduced for faster testing
        show_progress = TRUE,
        plot_vip = TRUE
    )
    
    cat("\n=== Analysis Results Summary ===\n")
    
    # Display omnibus results
    if (!is.null(results$omnibus)) {
        cat("\n--- Omnibus Analysis ---\n")
        cat("F-statistic:", round(results$omnibus$F_observed, 4), "\n")
        cat("P-value:", round(results$omnibus$p_value, 4), "\n")
        cat("Significant overall difference:", 
            ifelse(results$omnibus$p_value < 0.05, "YES", "NO"), "\n")
        
        # Show top globally important features
        top_global <- head(sort(results$omnibus$feature_weights, decreasing = TRUE), 10)
        cat("\nTop 10 globally important taxa:\n")
        for (i in 1:length(top_global)) {
            cat(sprintf("  %d. %s: %.4f\n", i, names(top_global)[i], top_global[i]))
        }
    }
    
    # Display pairwise results
    if (!is.null(results$pairwise)) {
        cat("\n--- Pairwise Analysis ---\n")
        cat("Total comparisons:", nrow(results$pairwise$summary_table), "\n")
        cat("Significant pairs (FDR corrected):", nrow(results$pairwise$significant_pairs), "\n")
        
        cat("\nAll pairwise comparisons:\n")
        print(results$pairwise$summary_table)
        
        if (nrow(results$pairwise$significant_pairs) > 0) {
            cat("\nSignificant pairwise comparisons:\n")
            for (i in 1:nrow(results$pairwise$significant_pairs)) {
                row <- results$pairwise$significant_pairs[i, ]
                cat(sprintf("  %s vs %s: F=%.3f, p=%.4f (corrected: %.4f)\n",
                           row$Group1, row$Group2, row$F_statistic, 
                           row$P_value, row$P_value_corrected))
            }
        }
        
        # Show feature importance for each significant pair
        if (nrow(results$pairwise$significant_pairs) > 0) {
            cat("\n--- Feature Importance by Comparison ---\n")
            for (i in 1:nrow(results$pairwise$significant_pairs)) {
                row <- results$pairwise$significant_pairs[i, ]
                comparison_name <- paste(row$Group1, row$Group2, sep = "_vs_")
                
                if (comparison_name %in% names(results$pairwise$pairwise_results)) {
                    pair_result <- results$pairwise$pairwise_results[[comparison_name]]
                    top_features <- head(sort(pair_result$feature_weights, decreasing = TRUE), 5)
                    
                    cat(sprintf("\n%s vs %s - Top 5 important taxa:\n", row$Group1, row$Group2))
                    for (j in 1:length(top_features)) {
                        cat(sprintf("  %d. %s: %.4f\n", j, names(top_features)[j], top_features[j]))
                    }
                }
            }
        }
    }
    
    # Compare with ground truth
    cat("\n=== Ground Truth Comparison ===\n")
    cat("Expected signals:\n")
    cat("- Global progressive signal: Taxa", paste(data$signal_info$global_taxa, collapse = ", "), "\n")
    cat("- Control vs Others signal: Taxa", paste(data$signal_info$control_vs_others_taxa, collapse = ", "), "\n")
    if (length(data$signal_info$mild_vs_severe_taxa) > 0) {
        cat("- Mild vs Severe signal: Taxa", paste(data$signal_info$mild_vs_severe_taxa, collapse = ", "), "\n")
    }
    
    return(results)
}

# Function to test different signal strengths
test_signal_strengths <- function() {
    
    signal_levels <- c("weak", "medium", "strong")
    results_list <- list()
    
    for (strength in signal_levels) {
        cat("\n" %R% rep("=", 60) %R% "\n")
        cat("Testing signal strength:", strength, "\n")
        cat(rep("=", 60), "\n")
        
        results_list[[strength]] <- run_multi_group_test(signal_strength = strength)
        
        # Brief pause between tests
        Sys.sleep(2)
    }
    
    # Summary comparison
    cat("\n" %R% rep("=", 60) %R% "\n")
    cat("SUMMARY COMPARISON ACROSS SIGNAL STRENGTHS\n")
    cat(rep("=", 60), "\n")
    
    summary_table <- data.frame(
        Signal_Strength = signal_levels,
        Omnibus_F = numeric(length(signal_levels)),
        Omnibus_P = numeric(length(signal_levels)),
        Significant_Pairs = numeric(length(signal_levels)),
        stringsAsFactors = FALSE
    )
    
    for (i in 1:length(signal_levels)) {
        strength <- signal_levels[i]
        if (!is.null(results_list[[strength]]$omnibus)) {
            summary_table$Omnibus_F[i] <- results_list[[strength]]$omnibus$F_observed
            summary_table$Omnibus_P[i] <- results_list[[strength]]$omnibus$p_value
        }
        if (!is.null(results_list[[strength]]$pairwise)) {
            summary_table$Significant_Pairs[i] <- nrow(results_list[[strength]]$pairwise$significant_pairs)
        }
    }
    
    print(summary_table)
    
    return(results_list)
}

# Run the tests
if (interactive()) {
    cat("MeLSI Multi-Group Testing Suite\n")
    cat("==============================\n\n")
    
    # Test with medium signal strength first
    cat("Running single test with medium signal strength...\n")
    medium_results <- run_multi_group_test(signal_strength = "medium")
    
    # Ask if user wants to run full comparison
    cat("\nWould you like to run tests across all signal strengths? (y/n): ")
    run_full <- readline()
    
    if (tolower(run_full) %in% c("y", "yes")) {
        all_results <- test_signal_strengths()
        cat("\nAll tests completed! Results stored in 'all_results' object.\n")
    }
} else {
    # Non-interactive mode - just run medium signal test
    cat("Running non-interactive multi-group test...\n")
    medium_results <- run_multi_group_test(signal_strength = "medium")
}
