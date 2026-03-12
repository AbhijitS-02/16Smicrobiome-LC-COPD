#!/usr/bin/env Rscript

# optimize_dada2_params.R
# This script performs a grid search for DADA2 truncation and maxEE parameters.
# It uses the trimmed reads as input.
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (requireNamespace("dada2", quietly = TRUE)) {
  cat("dada2 is installed.\n")
} else {
  cat("Installing dada2.\n")
  BiocManager::install("dada2", ask = FALSE, update = FALSE)
}

library(dada2)
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)

# --- CONFIGURATION ---
# Adjust these paths as needed. 
# Using hardcoded paths based on the project structure.
project_dir <- "D:/LC_COPD_microbiome"
trimmed_dir <- file.path(project_dir, "results/trimmed_reads")
output_dir <- file.path(project_dir, "results/dada2_optimization")

# File patterns
r1_pattern <- "_R1_trimmed.fastq.gz"
r2_pattern <- "_R2_trimmed.fastq.gz"

# Create output directory
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# --- 1. GET FILE LIST ---
cat("Listing files from:", trimmed_dir, "\n")
fnFs <- sort(list.files(trimmed_dir, pattern = r1_pattern, full.names = TRUE))
fnRs <- sort(list.files(trimmed_dir, pattern = r2_pattern, full.names = TRUE))

# Extract sample names
sample.names <- sapply(strsplit(basename(fnFs), "_R1"), `[`, 1)

if (length(fnFs) == 0) {
  stop("No forward read files found.")
}
if (length(fnFs) != length(fnRs)) {
  stop("Mismatch between number of forward and reverse files.")
}

cat("Found", length(fnFs), "samples.\n")

# --- 2. DEFINE GRID ---
# Truncation lengths to test (around current 278/240)
# Testing slightly wider range to find optimum
trunc_f_opts <- c(260, 270, 275) 
trunc_r_opts <- c(210, 220)

# MaxEE options
max_ee_opts <- list(c(2, 4), c(3, 4)) # c(maxEE_F, maxEE_R)

# Create grid of all unique combinations
grid_params <- expand.grid(
  truncLen_F = trunc_f_opts,
  truncLen_R = trunc_r_opts,
  maxEE_idx = 1:length(max_ee_opts)
)

cat("Testing", nrow(grid_params), "parameter combinations.\n")

# --- 3. RUN GRID SEARCH ---
# Running on all samples as requested.
cat("Using ALL", length(fnFs), "samples for optimization.\n")

results_list <- list()

for (i in 1:nrow(grid_params)) {
  t_f <- grid_params$truncLen_F[i]
  t_r <- grid_params$truncLen_R[i]
  ee_idx <- grid_params$maxEE_idx[i]
  ee_f <- max_ee_opts[[ee_idx]][1]
  ee_r <- max_ee_opts[[ee_idx]][2]
  
  cat(sprintf("\nTesting combination %d/%d: TruncF=%d, TruncR=%d, EE=(%d,%d)\n", 
              i, nrow(grid_params), t_f, t_r, ee_f, ee_r))
  
  # Temporary filenames for filtered reads

  filtFs_tmp <- file.path(output_dir, "temp_filt", paste0(sample.names, "_F_filt.fastq.gz"))
  filtRs_tmp <- file.path(output_dir, "temp_filt", paste0(sample.names, "_R_filt.fastq.gz"))
  
  # Run filterAndTrim
  out <- filterAndTrim(fnFs, filtFs_tmp, fnRs, filtRs_tmp,
                       truncLen = c(t_f, t_r),
                       maxEE = c(ee_f, ee_r),
                       rm.phix = TRUE,
                       compress = TRUE,
                       multithread = FALSE, # Set to FALSE if running on Windows
                       verbose = FALSE)
  
  # Calculate stats
  reads_in <- sum(out[, "reads.in"])
  reads_out <- sum(out[, "reads.out"])
  pct_retained <- (reads_out / reads_in) * 100
  
  cat(sprintf("  Reads In: %d, Reads Out: %d, Retained: %.2f%%\n", 
              reads_in, reads_out, pct_retained))
  
  results_list[[i]] <- data.frame(
    Combination = i,
    truncLen_F = t_f,
    truncLen_R = t_r,
    maxEE_F = ee_f,
    maxEE_R = ee_r,
    Reads_In = reads_in,
    Reads_Out = reads_out,
    Pct_Retained = pct_retained
  )
  
  # Clean up temp files
  unlink(file.path(output_dir, "temp_filt"), recursive = TRUE)
}

results_df <- do.call(rbind, results_list)

# Save optimization results
write.csv(results_df, file.path(output_dir, "dada2_optimization_results.csv"), row.names = FALSE)

# --- 4. SELECT TOP 3 ---
# Rank by percentage retained
top3 <- results_df %>%
  arrange(desc(Pct_Retained)) %>%
  slice(1:3)

cat("\nTop 3 Parameter Combinations:\n")
print(top3)

# --- 5. LEARN ERRORS & PLOT ---
pdf(file.path(output_dir, "dada2_error_plots.pdf"), width = 10, height = 10)

for (k in 1:nrow(top3)) {
  # Get params
  t_f <- top3$truncLen_F[k]
  t_r <- top3$truncLen_R[k]
  ee_f <- top3$maxEE_F[k]
  ee_r <- top3$maxEE_R[k]
  combo_id <- top3$Combination[k]
  
  cat(sprintf("\nRunning learnErrors for Top #%d (Combo %d)...\n", k, combo_id))
  
  # Need to re-run filterAndTrim (standard approach is to filter first, then learn errors)
  # We will use the same subset of samples to save time, or we could use more.
  # Let's stick to the subset for the error learning visualization as requested "learn error rates for those top 3".
  
  filt_path <- file.path(output_dir, paste0("filtered_top", k))
  if(!dir.exists(filt_path)) dir.create(filt_path)
  
  filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
  filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))
  
  out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                       truncLen = c(t_f, t_r),
                       maxEE = c(ee_f, ee_r),
                       rm.phix = TRUE,
                       compress = TRUE,
                       multithread = FALSE, # Set to FALSE if running on Windows
                       verbose = FALSE)
  
  # Learn errors
  errF <- learnErrors(filtFs, multithread=TRUE, verbose=FALSE)
  errR <- learnErrors(filtRs, multithread=TRUE, verbose=FALSE)
  
  # Plot
  p1 <- plotErrors(errF, nominalQ=TRUE) + 
    ggtitle(sprintf("Forward Errors (Top #%d: Trunc=%d, EE=%d)", k, t_f, ee_f))
  
  p2 <- plotErrors(errR, nominalQ=TRUE) + 
    ggtitle(sprintf("Reverse Errors (Top #%d: Trunc=%d, EE=%d)", k, t_r, ee_r))
  
  print(grid.arrange(p1, p2, nrow=2))
  
  cat("  Plots generated.\n")
  
  # Clean up filtered files to save space? 
  # Maybe keep them if the user wants to inspect, but for now we deletes to keep clean.
  # unlink(filt_path, recursive = TRUE)
}

dev.off()
cat("\nError plots saved to:", file.path(output_dir, "dada2_error_plots.pdf"), "\n")
cat("Optimization complete.\n")
