#' Generate Synthetic Microbiome Data for Testing
#'
#' Creates synthetic microbiome count data with specified characteristics for
#' testing and demonstration purposes.
#'
#' @param n_samples Total number of samples to generate
#' @param n_taxa Total number of taxa (features)
#' @param n_signal_taxa Number of taxa with differential abundance between groups
#' @param effect_size Magnitude of differential effect (fold-change)
#' @param group_balance Proportion of samples in group A (default 0.5 for balanced)
#'
#' @return A list containing:
#' \itemize{
#'   \item counts: Matrix of count data (samples x taxa)
#'   \item metadata: Data frame with sample metadata including Group
#'   \item signal_taxa: Vector of indices for taxa with true differential abundance
#' }
#'
#' @export
#'
#' @examples
#' # Generate balanced dataset with 60 samples and 100 taxa
#' test_data <- generate_test_data(n_samples = 60, n_taxa = 100, n_signal_taxa = 10)
#' X <- test_data$counts
#' y <- test_data$metadata$Group
#'
generate_test_data <- function(n_samples = 60,
                               n_taxa = 100,
                               n_signal_taxa = 10,
                               effect_size = 2.0,
                               group_balance = 0.5) {
  
  # Create group labels
  n_group_a <- round(n_samples * group_balance)
  n_group_b <- n_samples - n_group_a
  groups <- c(rep("A", n_group_a), rep("B", n_group_b))
  
  # Generate base count matrix (Poisson-distributed)
  base_lambda <- 100
  counts <- matrix(
    rpois(n_samples * n_taxa, lambda = base_lambda),
    nrow = n_samples,
    ncol = n_taxa
  )
  
  # Add differential abundance to signal taxa
  signal_taxa_idx <- seq_len(n_signal_taxa)
  
  for (i in signal_taxa_idx) {
    # Increase abundance in group B for signal taxa
    group_b_idx <- which(groups == "B")
    counts[group_b_idx, i] <- rpois(
      length(group_b_idx),
      lambda = base_lambda * effect_size
    )
  }
  
  # Add some sparsity (zeros) randomly
  zero_prop <- 0.2
  n_zeros <- round(length(counts) * zero_prop)
  zero_idx <- sample(length(counts), n_zeros)
  counts[zero_idx] <- 0
  
  # Create sample and taxa names
  rownames(counts) <- paste0("Sample", seq_len(n_samples))
  colnames(counts) <- paste0("Taxa", seq_len(n_taxa))
  
  # Create metadata
  metadata <- data.frame(
    SampleID = rownames(counts),
    Group = groups,
    stringsAsFactors = FALSE
  )
  rownames(metadata) <- metadata$SampleID
  
  # Return results
  return(list(
    counts = counts,
    metadata = metadata,
    signal_taxa = signal_taxa_idx
  ))
}




