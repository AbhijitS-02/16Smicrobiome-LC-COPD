# ===============================================================================
# 08_lefse.R - Linear discriminant analysis Effect Size (LEfSe)
# ===============================================================================
# Identify biomarker taxa using LEfSe
# Comparisons: Control, LC_COPD, LC_woCOPD
# ===============================================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(ggplot2)
  library(dplyr)
  # library(microbiomeMarker) # Required for LEfSe
})

# NOTE: If microbiomeMarker is not installed, install it using:
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install("microbiomeMarker")
library(microbiomeMarker)

# Using the project's custom plotting theme if available
source_path <- "D:/LC_COPD_microbiome/R_analysis/functions/plotting_theme.R"
if (file.exists(source_path)) {
  source(source_path)
}

# Configuration
project_dir <- "D:/LC_COPD_microbiome"
data_dir <- file.path(project_dir, "R_analysis/data")
output_dir <- file.path(project_dir, "figures/lefse")
tables_dir <- file.path(project_dir, "results/tables/lefse")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

cat("============================================\n")
cat("LEfSe Biomarker Discovery Analysis\n")
cat("============================================\n\n")

# Load data
cat("Loading phyloseq object...\n")
ps <- readRDS(file.path(data_dir, "phyloseq_silva_filtered.rds"))

# LEfSe is computationally intensive. We often perform it at a specific taxonomic
# rank (e.g., Genus) to find clinically relevant biomarkers.
cat("Aggregating taxa to Genus level...\n")
suppressMessages({
  ps_genus <- tax_glom(ps, taxrank = "Genus", NArm = FALSE)
})
table(sample_data(ps_genus)$group)

# Run LEfSe
cat("Running LEfSe analysis. This may take a few minutes...\n")
# LEfSe compares based on a class (our 'group' column)
set.seed(42)

# microbiomeMarker provides the run_lefse function:
# We specify the phyloseq object and the class variable.
# We use standard typical thresholds (wilcoxon p-value < 0.05, kw p-value < 0.05, lda score > 1)
lefse_out <- run_lefse(
  ps_genus,
  wilcoxon_cutoff = 0.05,
  norm = "CPM", # Common normalization for LEfSe (Counts per million)
  group = "group",
  kw_cutoff = 0.05,
  multigrp_strat = TRUE,
  lda_cutoff = 1.0
)

# Extract marker table
cat("Extracting results and saving tables...\n")
marker_df <- tryCatch(
  {
    as.data.frame(marker_table(lefse_out))
  },
  error = function(e) {
    message("No significant biomarkers detected.")
    return(data.frame())
  }
)

# Save result to CSV
write.csv(marker_df, file.path(tables_dir, "lefse_results_genus.csv"),
  row.names = FALSE
)

cat(sprintf("Found %d significant biomarkers (LDA > 1.0).\n", nrow(marker_df)))

# Generate Plotting
if (nrow(marker_df) > 0) {
  cat("Generating LDA bar plot...\n")

  # microbiomeMarker provides plot_ef_bar()
  p_bar <- plot_ef_bar(lefse_out) +
    theme_bw() +
    theme(
      axis.text.y = element_text(size = 8),
      legend.position = "top"
    ) +
    labs(title = "LEfSe LDA Scores", y = "LDA Score (log 10)")

  # Save bar plot
  ggsave(file.path(output_dir, "lefse_lda_barplot.png"), p_bar, width = 10, height = max(6, nrow(marker_df) * 0.15), dpi = 300)
  ggsave(file.path(output_dir, "lefse_lda_barplot.pdf"), p_bar, width = 10, height = max(6, nrow(marker_df) * 0.15))

  # Generate cladogram
  cat("Generating Cladogram...\n")
  # color vector should match the groups if specified, here we use generic
  p_clado <- as(plot_cladogram(lefse_out, color = c("Control" = "#195687ff", "LC_COPD" = "#ab3626ff", "LC_woCOPD" = "#00a087ff")), "ggplot")

  ggsave(file.path(output_dir, "lefse_cladogram.png"), p_clado, width = 12, height = 10, dpi = 300, bg = "white")
  ggsave(file.path(output_dir, "lefse_cladogram.pdf"), p_clado, width = 12, height = 10, bg = "white")
} else {
  cat("No significant biomarkers found to plot.\n")
}

cat("LEfSe analysis complete!\n")
