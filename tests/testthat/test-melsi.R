# Test file for MeLSI package
# This file contains basic tests for MeLSI functionality

library(testthat)
library(vegan)

# Source MeLSI functions
source("../../R/melsi_robust.R")

test_that("MeLSI basic functionality works", {
  
  # Generate test data
  test_data <- generate_test_data(n_samples = 40, n_taxa = 50, n_signal_taxa = 5)
  X <- test_data$counts
  y <- test_data$metadata$Group
  
  # CLR transformation
  X_clr <- X
  X_clr[X_clr == 0] <- 1e-10
  X_clr <- log(X_clr)
  X_clr <- X_clr - rowMeans(X_clr)
  
  # Test basic MeLSI run
  expect_no_error({
    results <- run_melsi_permutation_test(
      X_clr, y, 
      n_perms = 19, 
      B = 10, 
      show_progress = FALSE
    )
  })
  
  # Test results structure
  expect_true(is.list(results))
  expect_true("F_observed" %in% names(results))
  expect_true("p_value" %in% names(results))
  expect_true("M_ensemble" %in% names(results))
  
  # Test F-statistic is positive
  expect_true(results$F_observed > 0)
  
  # Test p-value is between 0 and 1
  expect_true(results$p_value >= 0 && results$p_value <= 1)
  
  # Test metric matrix is square
  expect_equal(nrow(results$M_ensemble), ncol(results$M_ensemble))
  expect_equal(ncol(results$M_ensemble), ncol(X_clr))
})

test_that("MeLSI outperforms Euclidean distance", {
  
  # Generate test data with known signal
  test_data <- generate_test_data(n_samples = 50, n_taxa = 30, n_signal_taxa = 6)
  X <- test_data$counts
  y <- test_data$metadata$Group
  
  # CLR transformation
  X_clr <- X
  X_clr[X_clr == 0] <- 1e-10
  X_clr <- log(X_clr)
  X_clr <- X_clr - rowMeans(X_clr)
  
  # Run MeLSI
  melsi_results <- run_melsi_permutation_test(
    X_clr, y, 
    n_perms = 19, 
    B = 20, 
    show_progress = FALSE
  )
  
  # Run Euclidean baseline
  dist_euclidean <- dist(X_clr)
  euclidean_results <- adonis2(dist_euclidean ~ y, permutations = 19)
  
  # MeLSI should have higher F-statistic
  expect_true(melsi_results$F_observed >= euclidean_results$F[1])
})

test_that("MeLSI handles edge cases", {
  
  # Test with minimal data
  X_minimal <- matrix(c(1, 2, 3, 4), nrow = 2, ncol = 2)
  y_minimal <- factor(c("A", "B"))
  
  # CLR transformation
  X_minimal_clr <- X_minimal
  X_minimal_clr[X_minimal_clr == 0] <- 1e-10
  X_minimal_clr <- log(X_minimal_clr)
  X_minimal_clr <- X_minimal_clr - rowMeans(X_minimal_clr)
  
  expect_no_error({
    results <- run_melsi_permutation_test(
      X_minimal_clr, y_minimal, 
      n_perms = 9, 
      B = 5, 
      show_progress = FALSE
    )
  })
  
  expect_true(results$F_observed > 0)
})
