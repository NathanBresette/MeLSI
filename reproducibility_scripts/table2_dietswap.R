#!/usr/bin/env Rscript
# ============================================================================
# Table 2 (DietSwap only): Method Comparison on DietSwap Dataset
# ============================================================================
# Reproduces the DietSwap portion of Table 2 using the real dataset
# provided in the microbiome R package. Compares MeLSI against
# Euclidean, Bray-Curtis, and Jaccard PERMANOVA analyses.
# ============================================================================

library(vegan)
library(MeLSI)
library(microbiome)

set.seed(42)

cat("\n==========================================================================\n")
cat("REPRODUCING TABLE 2 (DietSwap only)\n")
cat("==========================================================================\n\n")

# Load DietSwap dataset and subset
cat("--- Loading DietSwap dataset ---\n")
data(dietswap)
dietswap_subset <- subset_samples(dietswap, timepoint.within.group == 1)
dietswap_subset <- subset_samples(dietswap_subset, group %in% c("DI", "HE"))
dietswap_subset <- prune_taxa(taxa_sums(dietswap_subset) > 0, dietswap_subset)

counts_diet <- as(otu_table(dietswap_subset), "matrix")
if (taxa_are_rows(dietswap_subset)) {
  counts_diet <- t(counts_diet)
}
metadata_diet <- data.frame(sample_data(dietswap_subset))
metadata_diet$group <- droplevels(metadata_diet$group)
counts_diet <- counts_diet[rownames(metadata_diet), ]

cat(sprintf("Samples: %d\nTaxa: %d\nGroups: %s\n\n",
            nrow(counts_diet), ncol(counts_diet), paste(levels(metadata_diet$group), collapse = ", ")))

# CLR transform for MeLSI and Euclidean
X_diet_clr <- log(counts_diet + 1)
X_diet_clr <- X_diet_clr - rowMeans(X_diet_clr)
colnames(X_diet_clr) <- colnames(counts_diet)

y_diet <- metadata_diet$group

# Run MeLSI
cat("--- Running MeLSI on DietSwap ---\n")
melsi_diet <- melsi(X_diet_clr, y_diet, n_perms = 200, B = 30,
                    show_progress = TRUE, plot_vip = FALSE)

# Euclidean PERMANOVA
cat("--- Running Euclidean PERMANOVA ---\n")
dist_euc_diet <- dist(X_diet_clr)
perm_euc_diet <- adonis2(dist_euc_diet ~ y_diet, permutations = 999)

# Bray-Curtis PERMANOVA
cat("--- Running Bray-Curtis PERMANOVA ---\n")
dist_bray_diet <- vegdist(counts_diet, method = "bray")
perm_bray_diet <- adonis2(dist_bray_diet ~ y_diet, permutations = 999)

# Jaccard PERMANOVA
cat("--- Running Jaccard PERMANOVA ---\n")
dist_jaccard_diet <- vegdist(counts_diet, method = "jaccard", binary = TRUE)
perm_jaccard_diet <- adonis2(dist_jaccard_diet ~ y_diet, permutations = 999)

# DietSwap does not include a phylogenetic tree; UniFrac metrics unavailable
traditional_F <- c(
  Euclidean = perm_euc_diet$F[1],
  BrayCurtis = perm_bray_diet$F[1],
  Jaccard = perm_jaccard_diet$F[1]
)
best_trad <- names(which.max(traditional_F))

results <- data.frame(
  Dataset = "DietSwap (Real)",
  MeLSI_F = melsi_diet$F_observed,
  MeLSI_p = melsi_diet$p_value,
  Euclidean_F = perm_euc_diet$F[1],
  Euclidean_p = perm_euc_diet$`Pr(>F)`[1],
  BrayCurtis_F = perm_bray_diet$F[1],
  BrayCurtis_p = perm_bray_diet$`Pr(>F)`[1],
  Jaccard_F = perm_jaccard_diet$F[1],
  Jaccard_p = perm_jaccard_diet$`Pr(>F)`[1],
  WeightedUniFrac_F = NA,
  WeightedUniFrac_p = NA,
  UnweightedUniFrac_F = NA,
  UnweightedUniFrac_p = NA,
  Best_Traditional = best_trad,
  Best_Traditional_F = max(traditional_F)
)

print(results)
write.csv(results, "table2_dietswap_results.csv", row.names = FALSE)

cat("\nResults saved to table2_dietswap_results.csv\n")
cat("==========================================================================\n")
