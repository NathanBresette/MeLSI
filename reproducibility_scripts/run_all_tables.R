#!/usr/bin/env Rscript
# ==============================================================================
# Run All Tables for MeLSI Paper Reproducibility
# ==============================================================================
# This script runs all table reproduction scripts sequentially.
# ==============================================================================

cat("\n")
cat("================================================================================\n")
cat("MeLSI PAPER - COMPLETE REPRODUCIBILITY ANALYSIS\n")
cat("================================================================================\n")
cat("\n")
cat("This will reproduce all tables from the MeLSI paper.\n")
cat("Total estimated runtime: 60-90 minutes\n")
cat("\n")

start_time <- Sys.time()

# Table 1: Type I Error Control
cat("\n[1/5] Running Table 1: Type I Error Control\n")
cat("--------------------------------------------------------------------------------\n")
source("table1_type1_error.R")

# Table 2: Method Comparison on Synthetic and Real Datasets
cat("\n[2/5] Running Table 2: Method Comparison (Synthetic + Real)\n")
cat("--------------------------------------------------------------------------------\n")
source("table2_power_analysis.R")

# Table 3: Scalability
cat("\n[3/5] Running Table 3: Scalability\n")
cat("--------------------------------------------------------------------------------\n")
source("table3_scalability.R")

# Table 4: Parameter Sensitivity
cat("\n[4/5] Running Table 4: Parameter Sensitivity\n")
cat("--------------------------------------------------------------------------------\n")
source("table4_parameter_sensitivity.R")

# Table 5: Pre-filtering
cat("\n[5/5] Running Table 5: Pre-filtering Benefits\n")
cat("--------------------------------------------------------------------------------\n")
source("table5_prefiltering.R")

cat("\n")
cat("================================================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("================================================================================\n")
cat("Total runtime:", round(difftime(Sys.time(), start_time, units = "mins"), 1), "minutes\n")
cat("\nAll results saved to CSV files in this directory:\n")
cat("  - table1_results.csv\n")
cat("  - table2_results.csv\n")
cat("  - table3a_results.csv & table3b_results.csv\n")
cat("  - table4a_results.csv & table4b_results.csv\n")
cat("  - table5_results.csv\n")
cat("================================================================================\n")

