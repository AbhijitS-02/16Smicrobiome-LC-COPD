# ===============================================================================
# 05_ANCOMBC2.R - Differential Abundance Analysis with ANCOM-BC2
# ===============================================================================
# Identify differentially abundant taxa between groups
# Comparisons: LC_COPD vs Control, LC_woCOPD vs Control, LC_COPD vs LC_woCOPD
# Results are reported using two significance filters:
#   Section A: Raw p-value (p < 0.05)
#   Section B: Adjusted q-value / FDR (q < 0.05)
# Each filter gets its own output directory for figures and tables.
# Both combined and individual plots are exported.
# ===============================================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(ANCOMBC)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})

source("D:/LC_COPD_microbiome/R_analysis/functions/plotting_theme.R")

# Configuration
project_dir <- "D:/LC_COPD_microbiome"
data_dir <- file.path(project_dir, "R_analysis/data")

# Base output directories – sub-folders created later per filter
output_dir_base <- file.path(project_dir, "figures/differential_abundance")
tables_dir_base <- file.path(project_dir, "results/tables")

cat("============================================\n")
cat("ANCOM-BC2 Differential Abundance Analysis\n")
cat("============================================\n\n")

#-------------------------------------------------------------------------------
# LOAD AND PREPARE DATA
#-------------------------------------------------------------------------------
cat("Loading data...\n")
ps <- readRDS(file.path(data_dir, "phyloseq_silva_filtered.rds"))

# Filter low-abundance taxa for ANCOM-BC2
# Keep taxa present in at least 10% of samples
prevalence_threshold <- 0.10
prevalence <- apply(otu_table(ps), 1, function(x) sum(x > 0) / length(x))
ps_filtered <- prune_taxa(prevalence >= prevalence_threshold, ps)

cat(sprintf("  Taxa after filtering: %d (from %d)\n", ntaxa(ps_filtered), ntaxa(ps)))

# Aggregate to genus level for main analysis
ps_genus <- tax_glom(ps_filtered, taxrank = "Genus", NArm = FALSE)
cat(sprintf("  Genera for analysis: %d\n", ntaxa(ps_genus)))

#-------------------------------------------------------------------------------
# RUN ANCOM-BC2
#-------------------------------------------------------------------------------
cat("\nRunning ANCOM-BC2 (this may take a few minutes)...\n")

# Set reference group
sample_data(ps_genus)$group <- factor(
  sample_data(ps_genus)$group,
  levels = c("Control", "LC_COPD", "LC_woCOPD")
)

# Run ANCOM-BC2
set.seed(42)
ancom_output <- ancombc2(
  data = ps_genus,
  fix_formula = "group",
  p_adj_method = "BH",
  prv_cut = 0.10,
  lib_cut = 1000,
  group = "group",
  struc_zero = TRUE,
  neg_lb = TRUE,
  alpha = 0.05,
  global = TRUE,
  pairwise = TRUE,
  dunnet = TRUE,
  trend = FALSE,
  iter_control = list(tol = 1e-5, max_iter = 20, verbose = FALSE),
  em_control = list(tol = 1e-5, max_iter = 100),
  mdfdr_control = list(fwer_ctrl_method = "BH", B = 100),
  n_cl = 1
)

cat("  ANCOM-BC2 completed.\n")

#-------------------------------------------------------------------------------
# EXTRACT RESULTS
#-------------------------------------------------------------------------------
cat("\nExtracting results...\n")

# Get taxonomy for labeling
tax_df <- as.data.frame(tax_table(ps_genus))
tax_df$taxon <- rownames(tax_df)
tax_df$Taxon_Label <- ifelse(is.na(tax_df$Genus),
  paste0("Unclassified_", tax_df$Family),
  tax_df$Genus
)

# Primary results (vs Control)
res_primary <- ancom_output$res %>%
  left_join(tax_df %>% select(taxon, Taxon_Label, Phylum, Family, Genus),
    by = "taxon"
  )

# Pairwise results
res_pairwise <- ancom_output$res_pair %>%
  left_join(tax_df %>% select(taxon, Taxon_Label, Phylum, Family, Genus),
    by = "taxon"
  )

cat("  Primary result columns:\n")
cat("   ", paste(colnames(res_primary), collapse = ", "), "\n")
cat("  Pairwise result columns:\n")
cat("   ", paste(colnames(res_pairwise), collapse = ", "), "\n")

# ===============================================================================
# HELPER FUNCTIONS (shared by both filter sections)
# ===============================================================================

#--- Save individual + combined helper -----------------------------------------
save_plot <- function(plot_obj, dir, filename, width, height) {
  ggsave(file.path(dir, paste0(filename, ".pdf")),
    plot_obj,
    width = width, height = height, dpi = 300
  )
  ggsave(file.path(dir, paste0(filename, ".png")),
    plot_obj,
    width = width, height = height, dpi = 300
  )
}

#--- Volcano plot ---------------------------------------------------------------
create_volcano <- function(data, lfc_col, pval_col, sig_col_name, title,
                           yaxis_label = "-log10(value)") {
  df <- data %>%
    rename(
      lfc = !!sym(lfc_col),
      pval = !!sym(pval_col),
      diff = !!sym(sig_col_name)
    ) %>%
    mutate(
      Significance = case_when(
        diff & lfc > 0 ~ "Enriched",
        diff & lfc < 0 ~ "Depleted",
        TRUE ~ "Not significant"
      ),
      neg_log_p = -log10(pval)
    )

  # Label top taxa symmetrically — top 8 enriched + top 8 depleted
  top_enriched <- df %>%
    filter(diff & lfc > 0) %>%
    arrange(pval) %>%
    slice_head(n = 8)
  top_depleted <- df %>%
    filter(diff & lfc < 0) %>%
    arrange(pval) %>%
    slice_head(n = 8)
  top_taxa <- bind_rows(top_enriched, top_depleted)

  p <- ggplot(df, aes(x = lfc, y = neg_log_p, color = Significance)) +
    geom_point(alpha = 0.7, size = 2) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray40") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
    scale_color_manual(values = c(
      "Enriched" = "#ab3626ff",
      "Depleted" = "#195687ff",
      "Not significant" = "gray70"
    )) +
    ggrepel::geom_text_repel(
      data = top_taxa,
      aes(label = Taxon_Label),
      size = 3,
      max.overlaps = 15,
      box.padding = 0.5
    ) +
    labs(
      x = "Log2 Fold Change",
      y = yaxis_label,
      title = title
    ) +
    theme_publication() +
    theme(legend.position = "right")

  return(p)
}

#--- LFC bar plot (with significance stars) ------------------------------------
create_lfc_barplot <- function(data, lfc_col, pval_col, sig_col_name,
                               title, n_top = 20) {
  df <- data %>%
    rename(
      lfc = !!sym(lfc_col),
      pval = !!sym(pval_col),
      diff = !!sym(sig_col_name)
    ) %>%
    filter(diff) %>%
    arrange(desc(abs(lfc))) %>%
    slice_head(n = n_top) %>%
    mutate(
      Direction = ifelse(lfc > 0, "Enriched", "Depleted"),
      sig_label = case_when(
        pval < 0.001 ~ "***",
        pval < 0.01 ~ "**",
        pval < 0.05 ~ "*",
        TRUE ~ ""
      ),
      star_x = ifelse(lfc > 0, lfc + max(abs(lfc)) * 0.03,
        lfc - max(abs(lfc)) * 0.03
      ),
      star_hjust = ifelse(lfc > 0, 0, 1)
    ) %>%
    mutate(
      Taxon_Label = make.unique(Taxon_Label, sep = " "),
      Taxon_Label = factor(Taxon_Label, levels = Taxon_Label[order(lfc)])
    )

  if (nrow(df) == 0) {
    return(ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "No significant taxa") +
      theme_void())
  }

  p <- ggplot(df, aes(x = lfc, y = Taxon_Label, fill = Direction)) +
    geom_col(alpha = 0.85) +
    geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
    geom_text(aes(x = star_x, label = sig_label),
      size = 4, fontface = "bold", hjust = df$star_hjust
    ) +
    scale_fill_manual(values = c("Enriched" = "#ab3626ff", "Depleted" = "#195687ff")) +
    labs(x = "Log2 Fold Change", y = "", title = title) +
    theme_publication() +
    theme(
      axis.text.y = element_text(size = 9),
      legend.position = "top"
    )

  return(p)
}

#--- Waterfall plot -------------------------------------------------------------
create_waterfall_plot <- function(data, lfc_col, pval_col, sig_col_name,
                                  title, enriched_label = "Enriched",
                                  depleted_label = "Depleted") {
  df <- data %>%
    rename(
      lfc = !!sym(lfc_col),
      pval = !!sym(pval_col),
      diff = !!sym(sig_col_name)
    ) %>%
    filter(!is.na(lfc) & is.finite(lfc))

  # Keep top 15 enriched + top 15 depleted by LFC magnitude
  top_up <- df %>%
    filter(lfc > 0) %>%
    arrange(desc(lfc)) %>%
    slice_head(n = 15)
  top_down <- df %>%
    filter(lfc < 0) %>%
    arrange(lfc) %>%
    slice_head(n = 15)
  df <- bind_rows(top_up, top_down) %>%
    mutate(
      Direction = ifelse(lfc > 0, enriched_label, depleted_label),
      sig_label = case_when(
        pval < 0.001 ~ "***",
        pval < 0.01 ~ "**",
        pval < 0.05 ~ "*",
        TRUE ~ ""
      ),
      star_x = ifelse(lfc > 0, lfc + max(abs(lfc)) * 0.03,
        lfc - max(abs(lfc)) * 0.03
      ),
      star_hjust = ifelse(lfc > 0, 0, 1)
    ) %>%
    arrange(lfc) %>%
    mutate(
      Taxon_Label = make.unique(Taxon_Label, sep = " "),
      Taxon_Label = factor(Taxon_Label, levels = Taxon_Label)
    )

  if (nrow(df) == 0) {
    return(ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "No taxa to plot") +
      theme_void())
  }

  p <- ggplot(df, aes(x = lfc, y = Taxon_Label, fill = Direction)) +
    geom_col(width = 0.75) +
    geom_vline(xintercept = 0, color = "black", linewidth = 0.4) +
    geom_text(aes(x = star_x, label = sig_label),
      size = 4, fontface = "bold", hjust = df$star_hjust
    ) +
    scale_fill_manual(
      values = setNames(c("#ab3626ff", "#195687ff"), c(enriched_label, depleted_label)),
      name = ""
    ) +
    labs(x = "Log Fold Change", y = "", title = title) +
    theme_publication() +
    theme(
      axis.text.y = element_text(size = 8),
      legend.position = "top",
      legend.text = element_text(size = 11),
      panel.grid.major.y = element_blank()
    )

  return(p)
}

#--- Heatmap helper -------------------------------------------------------------
create_heatmap <- function(sig_taxa, ps_genus_obj, tax_df_ref, res_primary_ref,
                           output_path, title_text) {
  if (length(sig_taxa) == 0) {
    cat("    No significant taxa — skipping heatmap.\n")
    return(invisible(NULL))
  }

  ps_sig <- prune_taxa(sig_taxa, ps_genus_obj)
  ps_sig_rel <- transform_sample_counts(ps_sig, function(x) x / sum(x) * 100)

  abund_mat <- as(otu_table(ps_sig_rel), "matrix")
  if (!taxa_are_rows(ps_sig_rel)) abund_mat <- t(abund_mat)

  rownames(abund_mat) <- make.unique(
    tax_df_ref[rownames(abund_mat), "Taxon_Label"],
    sep = "_"
  )

  log_mat <- log10(abund_mat + 0.01)
  log_mat[!is.finite(log_mat)] <- min(log_mat[is.finite(log_mat)], na.rm = TRUE)

  samp_data <- as.data.frame(sample_data(ps_sig))
  annotation_col <- data.frame(Group = samp_data$group, row.names = rownames(samp_data))

  do_cluster_rows <- nrow(log_mat) >= 2
  do_cluster_cols <- ncol(log_mat) >= 2

  tryCatch(
    {
      pdf(output_path,
        width = 12, height = max(6, length(sig_taxa) * 0.3 + 2)
      )
      pheatmap::pheatmap(
        log_mat,
        annotation_col = annotation_col,
        annotation_colors = list(Group = group_colors),
        color = viridis::viridis(100),
        cluster_rows = do_cluster_rows,
        cluster_cols = do_cluster_cols,
        fontsize_row = 8,
        fontsize_col = 7,
        main = title_text
      )
      dev.off()
      cat("    Heatmap saved.\n")
    },
    error = function(e) {
      cat(sprintf("    Warning: Heatmap generation failed (%s). Skipping.\n", e$message))
      try(dev.off(), silent = TRUE)
    }
  )
}

# ===============================================================================
# UTILITY: Add a significance flag column based on raw p-value < 0.05
# (ANCOM-BC2 only provides diff_* for q-value; we recreate for p-value)
# ===============================================================================
add_pval_diff_columns <- function(df) {
  # Detect all p_group* columns
  p_cols <- grep("^p_group", colnames(df), value = TRUE)
  for (pc in p_cols) {
    new_col <- sub("^p_", "diff_pval_", pc)
    df[[new_col]] <- !is.na(df[[pc]]) & df[[pc]] < 0.05
  }
  return(df)
}

res_primary <- add_pval_diff_columns(res_primary)
res_pairwise <- add_pval_diff_columns(res_pairwise)

# ---------------------------------------------------------------------------
# Dynamically detect pairwise column names for LC_COPD vs LC_woCOPD
# ANCOM-BC2 may name columns differently depending on the version / contrasts.
# ---------------------------------------------------------------------------
pair_lfc_col <- grep("^lfc_.*LC_COPD.*LC_woCOPD", colnames(res_pairwise), value = TRUE)[1]
if (is.na(pair_lfc_col)) {
  # Try reversed order
  pair_lfc_col <- grep("^lfc_.*LC_woCOPD.*LC_COPD", colnames(res_pairwise), value = TRUE)[1]
}

has_pairwise <- !is.na(pair_lfc_col)

if (has_pairwise) {
  # Derive matching p, q, se, diff column names from the lfc column
  pair_suffix <- sub("^lfc_", "", pair_lfc_col)
  pair_p_col <- paste0("p_", pair_suffix)
  pair_q_col <- paste0("q_", pair_suffix)
  pair_se_col <- paste0("se_", pair_suffix)
  pair_diff_col <- paste0("diff_", pair_suffix) # q-value flag
  pair_diff_pval_col <- paste0("diff_pval_", pair_suffix) # p-value flag
  cat(sprintf("  Detected pairwise columns (suffix: %s)\n", pair_suffix))
} else {
  cat("  WARNING: No pairwise LC_COPD vs LC_woCOPD columns found in res_pair.\n")
  cat("           Available columns: ", paste(colnames(res_pairwise), collapse = ", "), "\n")
}

# Pre-compute n_taxa for height calculations
n_taxa_copd <- sum(!is.na(res_primary$lfc_groupLC_COPD))
n_taxa_wocopd <- sum(!is.na(res_primary$lfc_groupLC_woCOPD))
if (has_pairwise) {
  n_taxa_pair <- sum(!is.na(res_pairwise[[pair_lfc_col]]))
}

# ╔═════════════════════════════════════════════════════════════════════════════╗
# ║  SECTION A: SIGNIFICANCE BASED ON RAW P-VALUE (p < 0.05)                  ║
# ╚═════════════════════════════════════════════════════════════════════════════╝
cat("\n\n========================================================\n")
cat("SECTION A: Results filtered by RAW P-VALUE (p < 0.05)\n")
cat("========================================================\n\n")

output_dir_pval <- file.path(output_dir_base, "pval_filtered")
tables_dir_pval <- file.path(tables_dir_base, "pval_filtered")
dir.create(output_dir_pval, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir_pval, recursive = TRUE, showWarnings = FALSE)

# Count significant taxa (p-value)
sig_lc_copd_p <- sum(res_primary$diff_pval_groupLC_COPD, na.rm = TRUE)
sig_lc_wocopd_p <- sum(res_primary$diff_pval_groupLC_woCOPD, na.rm = TRUE)

cat(sprintf("  Significant vs Control (p-value):\n"))
cat(sprintf("    LC_COPD vs Control: %d taxa\n", sig_lc_copd_p))
cat(sprintf("    LC_woCOPD vs Control: %d taxa\n", sig_lc_wocopd_p))

if (has_pairwise) {
  sig_pair_p <- sum(res_pairwise[[pair_diff_pval_col]], na.rm = TRUE)
  cat(sprintf("    LC_COPD vs LC_woCOPD: %d taxa\n", sig_pair_p))
}

#--- Volcano plots (p-value) ---------------------------------------------------
cat("\n  Creating volcano plots (p-value filter)...\n")

volcano_lc_copd_p <- create_volcano(
  res_primary,
  "lfc_groupLC_COPD", "p_groupLC_COPD", "diff_pval_groupLC_COPD",
  "LC_COPD vs Control",
  yaxis_label = "-log10(p-value)"
)

volcano_lc_wocopd_p <- create_volcano(
  res_primary,
  "lfc_groupLC_woCOPD", "p_groupLC_woCOPD", "diff_pval_groupLC_woCOPD",
  "LC_woCOPD vs Control",
  yaxis_label = "-log10(p-value)"
)

# Individual volcano plots
save_plot(volcano_lc_copd_p, output_dir_pval,
  "volcano_LC_COPD_vs_Control",
  width = 8, height = 6
)
save_plot(volcano_lc_wocopd_p, output_dir_pval,
  "volcano_LC_woCOPD_vs_Control",
  width = 8, height = 6
)

# Pairwise volcano (LC_COPD vs LC_woCOPD) — ALWAYS generated when columns exist
if (has_pairwise) {
  volcano_pair_p <- create_volcano(
    res_pairwise,
    pair_lfc_col,
    pair_p_col,
    pair_diff_pval_col,
    "LC_COPD vs LC_woCOPD",
    yaxis_label = "-log10(p-value)"
  )
  save_plot(volcano_pair_p, output_dir_pval,
    "volcano_LC_COPD_vs_LC_woCOPD",
    width = 8, height = 6
  )

  # Combined volcano (all three)
  volcano_combined_p <- volcano_lc_copd_p + volcano_lc_wocopd_p + volcano_pair_p +
    plot_layout(ncol = 3, guides = "collect") +
    plot_annotation(title = "Differential Abundance — ANCOM-BC2 (p-value filter)")
  save_plot(volcano_combined_p, output_dir_pval,
    "volcano_plots_combined",
    width = 20, height = 6
  )
} else {
  volcano_combined_p <- volcano_lc_copd_p + volcano_lc_wocopd_p +
    plot_layout(guides = "collect") +
    plot_annotation(title = "Differential Abundance — ANCOM-BC2 (p-value filter)")
  save_plot(volcano_combined_p, output_dir_pval,
    "volcano_plots_combined",
    width = 14, height = 6
  )
}

#--- LFC bar plots (p-value) ---------------------------------------------------
cat("  Creating bar plots (p-value filter)...\n")

bar_lc_copd_p <- create_lfc_barplot(
  res_primary,
  "lfc_groupLC_COPD", "p_groupLC_COPD", "diff_pval_groupLC_COPD",
  "LC_COPD vs Control"
)

bar_lc_wocopd_p <- create_lfc_barplot(
  res_primary,
  "lfc_groupLC_woCOPD", "p_groupLC_woCOPD", "diff_pval_groupLC_woCOPD",
  "LC_woCOPD vs Control"
)

# Individual bar plots
save_plot(bar_lc_copd_p, output_dir_pval,
  "lfc_barplot_LC_COPD_vs_Control",
  width = 8, height = 8
)
save_plot(bar_lc_wocopd_p, output_dir_pval,
  "lfc_barplot_LC_woCOPD_vs_Control",
  width = 8, height = 8
)

# Pairwise bar plot (LC_COPD vs LC_woCOPD)
if (has_pairwise) {
  bar_pair_p <- create_lfc_barplot(
    res_pairwise,
    pair_lfc_col,
    pair_p_col,
    pair_diff_pval_col,
    "LC_COPD vs LC_woCOPD"
  )
  save_plot(bar_pair_p, output_dir_pval,
    "lfc_barplot_LC_COPD_vs_LC_woCOPD",
    width = 8, height = 8
  )

  # Combined bar plots (all three)
  bar_combined_p <- bar_lc_copd_p + bar_lc_wocopd_p + bar_pair_p +
    plot_layout(ncol = 3, guides = "collect") &
    theme(legend.position = "top")
  save_plot(bar_combined_p, output_dir_pval,
    "lfc_barplots_combined",
    width = 20, height = 8
  )
} else {
  bar_combined_p <- bar_lc_copd_p + bar_lc_wocopd_p +
    plot_layout(guides = "collect") &
    theme(legend.position = "top")
  save_plot(bar_combined_p, output_dir_pval,
    "lfc_barplots_combined",
    width = 14, height = 8
  )
}

#--- Waterfall plots (p-value) -------------------------------------------------
cat("  Creating waterfall plots (p-value filter)...\n")

waterfall_copd_p <- create_waterfall_plot(
  res_primary,
  "lfc_groupLC_COPD", "p_groupLC_COPD", "diff_pval_groupLC_COPD",
  "LC_COPD vs Control",
  enriched_label = "Higher in LC_COPD",
  depleted_label = "Higher in Control"
)
save_plot(waterfall_copd_p, output_dir_pval,
  "waterfall_lfc_LC_COPD_vs_Control",
  width = 10, height = max(6, n_taxa_copd * 0.22)
)

waterfall_wocopd_p <- create_waterfall_plot(
  res_primary,
  "lfc_groupLC_woCOPD", "p_groupLC_woCOPD", "diff_pval_groupLC_woCOPD",
  "LC_woCOPD vs Control",
  enriched_label = "Higher in LC_woCOPD",
  depleted_label = "Higher in Control"
)
save_plot(waterfall_wocopd_p, output_dir_pval,
  "waterfall_lfc_LC_woCOPD_vs_Control",
  width = 10, height = max(6, n_taxa_wocopd * 0.22)
)

# Pairwise waterfall (p-value)
if (has_pairwise) {
  waterfall_pair_p <- create_waterfall_plot(
    res_pairwise,
    pair_lfc_col,
    pair_p_col,
    pair_diff_pval_col,
    "LC_COPD vs LC_woCOPD",
    enriched_label = "Higher in LC_COPD",
    depleted_label = "Higher in LC_woCOPD"
  )
  save_plot(waterfall_pair_p, output_dir_pval,
    "waterfall_lfc_LC_COPD_vs_LC_woCOPD",
    width = 10, height = max(6, n_taxa_pair * 0.22)
  )
}
cat("  Waterfall plots (p-value) saved.\n")

#--- Heatmap (p-value) ---------------------------------------------------------
cat("  Creating heatmap (p-value filter)...\n")

sig_taxa_pval <- res_primary %>%
  filter(diff_pval_groupLC_COPD | diff_pval_groupLC_woCOPD) %>%
  pull(taxon)

create_heatmap(
  sig_taxa_pval, ps_genus, tax_df, res_primary,
  file.path(output_dir_pval, "heatmap_significant_taxa.pdf"),
  "Differentially Abundant Taxa — ANCOM-BC2 (p-value < 0.05)"
)

#--- Summary tables (p-value) --------------------------------------------------
cat("  Creating summary tables (p-value filter)...\n")

summary_table_pval <- res_primary %>%
  select(
    taxon, Taxon_Label, Phylum, Family, Genus,
    lfc_groupLC_COPD, se_groupLC_COPD, p_groupLC_COPD, diff_pval_groupLC_COPD,
    lfc_groupLC_woCOPD, se_groupLC_woCOPD, p_groupLC_woCOPD, diff_pval_groupLC_woCOPD
  ) %>%
  arrange(p_groupLC_COPD, p_groupLC_woCOPD)

write.csv(summary_table_pval,
  file.path(tables_dir_pval, "ancombc2_results_pval.csv"),
  row.names = FALSE
)

sig_table_pval <- summary_table_pval %>%
  filter(diff_pval_groupLC_COPD | diff_pval_groupLC_woCOPD)

write.csv(sig_table_pval,
  file.path(tables_dir_pval, "ancombc2_significant_pval.csv"),
  row.names = FALSE
)

# Pairwise table (p-value)
if (has_pairwise) {
  summary_table_pair_pval <- res_pairwise %>%
    select(
      taxon, Taxon_Label, Phylum, Family, Genus,
      all_of(c(pair_lfc_col, pair_se_col, pair_p_col, pair_diff_pval_col))
    ) %>%
    arrange(.data[[pair_p_col]])

  write.csv(summary_table_pair_pval,
    file.path(tables_dir_pval, "ancombc2_results_pval_LC_COPD_vs_LC_woCOPD.csv"),
    row.names = FALSE
  )

  sig_table_pair_pval <- summary_table_pair_pval %>%
    filter(.data[[pair_diff_pval_col]])

  write.csv(sig_table_pair_pval,
    file.path(tables_dir_pval, "ancombc2_significant_pval_LC_COPD_vs_LC_woCOPD.csv"),
    row.names = FALSE
  )
}

cat(sprintf("  Tables saved (%d significant taxa by p-value).\n", nrow(sig_table_pval)))


# ╔═════════════════════════════════════════════════════════════════════════════╗
# ║  SECTION B: SIGNIFICANCE BASED ON Q-VALUE / FDR (q < 0.05)                ║
# ╚═════════════════════════════════════════════════════════════════════════════╝
cat("\n\n========================================================\n")
cat("SECTION B: Results filtered by ADJUSTED Q-VALUE (q < 0.05)\n")
cat("========================================================\n\n")

output_dir_qval <- file.path(output_dir_base, "qval_filtered")
tables_dir_qval <- file.path(tables_dir_base, "qval_filtered")
dir.create(output_dir_qval, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir_qval, recursive = TRUE, showWarnings = FALSE)

# Count significant taxa (q-value) — uses the original diff_* columns
sig_lc_copd_q <- sum(res_primary$diff_groupLC_COPD, na.rm = TRUE)
sig_lc_wocopd_q <- sum(res_primary$diff_groupLC_woCOPD, na.rm = TRUE)

cat(sprintf("  Significant vs Control (q-value):\n"))
cat(sprintf("    LC_COPD vs Control: %d taxa\n", sig_lc_copd_q))
cat(sprintf("    LC_woCOPD vs Control: %d taxa\n", sig_lc_wocopd_q))

if (has_pairwise) {
  sig_pair_q <- sum(res_pairwise[[pair_diff_col]], na.rm = TRUE)
  cat(sprintf("    LC_COPD vs LC_woCOPD: %d taxa\n", sig_pair_q))
}

#--- Volcano plots (q-value) ---------------------------------------------------
cat("\n  Creating volcano plots (q-value filter)...\n")

volcano_lc_copd_q <- create_volcano(
  res_primary,
  "lfc_groupLC_COPD", "q_groupLC_COPD", "diff_groupLC_COPD",
  "LC_COPD vs Control",
  yaxis_label = "-log10(q-value)"
)

volcano_lc_wocopd_q <- create_volcano(
  res_primary,
  "lfc_groupLC_woCOPD", "q_groupLC_woCOPD", "diff_groupLC_woCOPD",
  "LC_woCOPD vs Control",
  yaxis_label = "-log10(q-value)"
)

# Individual volcano plots
save_plot(volcano_lc_copd_q, output_dir_qval,
  "volcano_LC_COPD_vs_Control",
  width = 8, height = 6
)
save_plot(volcano_lc_wocopd_q, output_dir_qval,
  "volcano_LC_woCOPD_vs_Control",
  width = 8, height = 6
)

# Pairwise volcano (LC_COPD vs LC_woCOPD) — ALWAYS generated when columns exist
if (has_pairwise) {
  volcano_pair_q <- create_volcano(
    res_pairwise,
    pair_lfc_col,
    pair_q_col,
    pair_diff_col,
    "LC_COPD vs LC_woCOPD",
    yaxis_label = "-log10(q-value)"
  )
  save_plot(volcano_pair_q, output_dir_qval,
    "volcano_LC_COPD_vs_LC_woCOPD",
    width = 8, height = 6
  )

  # Combined volcano (all three)
  volcano_combined_q <- volcano_lc_copd_q + volcano_lc_wocopd_q + volcano_pair_q +
    plot_layout(ncol = 3, guides = "collect") +
    plot_annotation(title = "Differential Abundance — ANCOM-BC2 (q-value filter)")
  save_plot(volcano_combined_q, output_dir_qval,
    "volcano_plots_combined",
    width = 20, height = 6
  )
} else {
  volcano_combined_q <- volcano_lc_copd_q + volcano_lc_wocopd_q +
    plot_layout(guides = "collect") +
    plot_annotation(title = "Differential Abundance — ANCOM-BC2 (q-value filter)")
  save_plot(volcano_combined_q, output_dir_qval,
    "volcano_plots_combined",
    width = 14, height = 6
  )
}

#--- LFC bar plots (q-value) ---------------------------------------------------
cat("  Creating bar plots (q-value filter)...\n")

bar_lc_copd_q <- create_lfc_barplot(
  res_primary,
  "lfc_groupLC_COPD", "q_groupLC_COPD", "diff_groupLC_COPD",
  "LC_COPD vs Control"
)

bar_lc_wocopd_q <- create_lfc_barplot(
  res_primary,
  "lfc_groupLC_woCOPD", "q_groupLC_woCOPD", "diff_groupLC_woCOPD",
  "LC_woCOPD vs Control"
)

# Individual bar plots
save_plot(bar_lc_copd_q, output_dir_qval,
  "lfc_barplot_LC_COPD_vs_Control",
  width = 8, height = 8
)
save_plot(bar_lc_wocopd_q, output_dir_qval,
  "lfc_barplot_LC_woCOPD_vs_Control",
  width = 8, height = 8
)

# Pairwise bar plot (LC_COPD vs LC_woCOPD)
if (has_pairwise) {
  bar_pair_q <- create_lfc_barplot(
    res_pairwise,
    pair_lfc_col,
    pair_q_col,
    pair_diff_col,
    "LC_COPD vs LC_woCOPD"
  )
  save_plot(bar_pair_q, output_dir_qval,
    "lfc_barplot_LC_COPD_vs_LC_woCOPD",
    width = 8, height = 8
  )

  # Combined bar plots (all three)
  bar_combined_q <- bar_lc_copd_q + bar_lc_wocopd_q + bar_pair_q +
    plot_layout(ncol = 3, guides = "collect") &
    theme(legend.position = "top")
  save_plot(bar_combined_q, output_dir_qval,
    "lfc_barplots_combined",
    width = 20, height = 8
  )
} else {
  bar_combined_q <- bar_lc_copd_q + bar_lc_wocopd_q +
    plot_layout(guides = "collect") &
    theme(legend.position = "top")
  save_plot(bar_combined_q, output_dir_qval,
    "lfc_barplots_combined",
    width = 14, height = 8
  )
}

#--- Waterfall plots (q-value) -------------------------------------------------
cat("  Creating waterfall plots (q-value filter)...\n")

waterfall_copd_q <- create_waterfall_plot(
  res_primary,
  "lfc_groupLC_COPD", "q_groupLC_COPD", "diff_groupLC_COPD",
  "LC_COPD vs Control",
  enriched_label = "Higher in LC_COPD",
  depleted_label = "Higher in Control"
)
save_plot(waterfall_copd_q, output_dir_qval,
  "waterfall_lfc_LC_COPD_vs_Control",
  width = 10, height = max(6, n_taxa_copd * 0.22)
)

waterfall_wocopd_q <- create_waterfall_plot(
  res_primary,
  "lfc_groupLC_woCOPD", "q_groupLC_woCOPD", "diff_groupLC_woCOPD",
  "LC_woCOPD vs Control",
  enriched_label = "Higher in LC_woCOPD",
  depleted_label = "Higher in Control"
)
save_plot(waterfall_wocopd_q, output_dir_qval,
  "waterfall_lfc_LC_woCOPD_vs_Control",
  width = 10, height = max(6, n_taxa_wocopd * 0.22)
)

# Pairwise waterfall (q-value)
if (has_pairwise) {
  waterfall_pair_q <- create_waterfall_plot(
    res_pairwise,
    pair_lfc_col,
    pair_q_col,
    pair_diff_col,
    "LC_COPD vs LC_woCOPD",
    enriched_label = "Higher in LC_COPD",
    depleted_label = "Higher in LC_woCOPD"
  )
  save_plot(waterfall_pair_q, output_dir_qval,
    "waterfall_lfc_LC_COPD_vs_LC_woCOPD",
    width = 10, height = max(6, n_taxa_pair * 0.22)
  )
}
cat("  Waterfall plots (q-value) saved.\n")

#--- Heatmap (q-value) ---------------------------------------------------------
cat("  Creating heatmap (q-value filter)...\n")

sig_taxa_qval <- res_primary %>%
  filter(diff_groupLC_COPD | diff_groupLC_woCOPD) %>%
  pull(taxon)

create_heatmap(
  sig_taxa_qval, ps_genus, tax_df, res_primary,
  file.path(output_dir_qval, "heatmap_significant_taxa.pdf"),
  "Differentially Abundant Taxa — ANCOM-BC2 (q-value < 0.05)"
)

#--- Summary tables (q-value) --------------------------------------------------
cat("  Creating summary tables (q-value filter)...\n")

summary_table_qval <- res_primary %>%
  select(
    taxon, Taxon_Label, Phylum, Family, Genus,
    lfc_groupLC_COPD, se_groupLC_COPD, q_groupLC_COPD, diff_groupLC_COPD,
    lfc_groupLC_woCOPD, se_groupLC_woCOPD, q_groupLC_woCOPD, diff_groupLC_woCOPD
  ) %>%
  arrange(q_groupLC_COPD, q_groupLC_woCOPD)

write.csv(summary_table_qval,
  file.path(tables_dir_qval, "ancombc2_results_qval.csv"),
  row.names = FALSE
)

sig_table_qval <- summary_table_qval %>%
  filter(diff_groupLC_COPD | diff_groupLC_woCOPD)

write.csv(sig_table_qval,
  file.path(tables_dir_qval, "ancombc2_significant_qval.csv"),
  row.names = FALSE
)

# Pairwise table (q-value)
if (has_pairwise) {
  summary_table_pair_qval <- res_pairwise %>%
    select(
      taxon, Taxon_Label, Phylum, Family, Genus,
      all_of(c(pair_lfc_col, pair_se_col, pair_q_col, pair_diff_col))
    ) %>%
    arrange(.data[[pair_q_col]])

  write.csv(summary_table_pair_qval,
    file.path(tables_dir_qval, "ancombc2_results_qval_LC_COPD_vs_LC_woCOPD.csv"),
    row.names = FALSE
  )

  sig_table_pair_qval <- summary_table_pair_qval %>%
    filter(.data[[pair_diff_col]])

  write.csv(sig_table_pair_qval,
    file.path(tables_dir_qval, "ancombc2_significant_qval_LC_COPD_vs_LC_woCOPD.csv"),
    row.names = FALSE
  )
}

cat(sprintf("  Tables saved (%d significant taxa by q-value).\n", nrow(sig_table_qval)))


# ===============================================================================
# COMBINED SUMMARY
# ===============================================================================
cat("\n\n============================================\n")
cat("ANCOM-BC2 Results Summary\n")
cat("============================================\n")
cat(sprintf("Total taxa tested: %d\n", nrow(res_primary)))

cat("\n--- Raw p-value filter (p < 0.05) ---\n")
cat(sprintf("  Significant LC_COPD vs Control:   %d\n", sig_lc_copd_p))
cat(sprintf("  Significant LC_woCOPD vs Control:  %d\n", sig_lc_wocopd_p))
if (has_pairwise) {
  cat(sprintf("  Significant LC_COPD vs LC_woCOPD: %d\n", sig_pair_p))
}

cat("\n--- Adjusted q-value filter (q < 0.05) ---\n")
cat(sprintf("  Significant LC_COPD vs Control:   %d\n", sig_lc_copd_q))
cat(sprintf("  Significant LC_woCOPD vs Control:  %d\n", sig_lc_wocopd_q))
if (has_pairwise) {
  cat(sprintf("  Significant LC_COPD vs LC_woCOPD: %d\n", sig_pair_q))
}

if (nrow(sig_table_pval) > 0) {
  cat("\nTop significant taxa (p-value):\n")
  print(head(sig_table_pval %>% select(Taxon_Label, lfc_groupLC_COPD, p_groupLC_COPD), 10))
}

if (nrow(sig_table_qval) > 0) {
  cat("\nTop significant taxa (q-value):\n")
  print(head(sig_table_qval %>% select(Taxon_Label, lfc_groupLC_COPD, q_groupLC_COPD), 10))
}

# Save ANCOM-BC2 object
saveRDS(ancom_output, file.path(data_dir, "ancombc2_output.rds"))

cat("\n============================================\n")
cat("ANCOM-BC2 Analysis Complete!\n")
cat("============================================\n")
