# ===============================================================================
# 07_PUBLICATION_PLOTS.R - Generate Publication-Quality Figures
# ===============================================================================
# Create final publication-ready figures combining all analyses
#
# Figures produced:
#   Figure 1  – Alpha Diversity Panel (4 metrics)
#   Figure 2  – Beta Diversity Ordination (Bray-Curtis + wUniFrac)
#   Figure 3  – Taxonomic Composition (phylum-level stacked bar)
#   Figure 4  – Differential Abundance (ANCOM-BC2 LFC barplot)
#   Figure 5  – Functional Pathway Analysis (PCA + significant barplot)
#   Figure S1 – Rarefaction Curves (copy from existing)
#   Figure S2 – Volcano Plots (combined ANCOM-BC2)
#   Figure S3 – Genus-Level Heatmap + Core Microbiome Venn
#   Figure S4 – KO & EC Significant Barplots
#   Table 1   – Sample Summary Statistics
#
# All PNG outputs saved at 600 DPI.
# ===============================================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(ggplot2)
  library(ggpubr)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(grid)
  library(gridExtra)
  library(vegan)
  library(pheatmap)
  library(viridis)
  library(ggrepel)
})

source("D:/LC_COPD_microbiome/R_analysis/functions/plotting_theme.R")

# Configuration
project_dir <- "D:/LC_COPD_microbiome"
data_dir <- file.path(project_dir, "R_analysis/data")
output_dir <- file.path(project_dir, "figures")
pub_dir <- file.path(output_dir, "publication")
tables_dir <- file.path(project_dir, "results/tables")
metadata_file <- file.path(project_dir, "metadata/sample_metadata.tsv")

dir.create(pub_dir, recursive = TRUE, showWarnings = FALSE)

# Global DPI for all PNG saves
PUB_DPI <- 600

cat("============================================\n")
cat("Generating Publication-Quality Figures\n")
cat("============================================\n\n")

# ---------------------------------------------------------------------------
# Helper: save plot in PDF + PNG (600 DPI) + TIFF
# ---------------------------------------------------------------------------
save_pub_plot <- function(plot_obj, basename, width, height,
                          formats = c("pdf", "png", "tiff")) {
  for (fmt in formats) {
    outpath <- file.path(pub_dir, paste0(basename, ".", fmt))
    if (fmt == "tiff") {
      ggsave(outpath, plot_obj,
        width = width, height = height, dpi = PUB_DPI,
        compression = "lzw"
      )
    } else {
      ggsave(outpath, plot_obj,
        width = width, height = height, dpi = PUB_DPI
      )
    }
  }
  cat(sprintf("  Saved: %s (.pdf/.png/.tiff)\n", basename))
}

#-------------------------------------------------------------------------------
# LOAD ALL DATA
#-------------------------------------------------------------------------------
cat("Loading data...\n")

# Phyloseq object
ps <- readRDS(file.path(data_dir, "phyloseq_silva_filtered.rds"))
ps_rel <- transform_sample_counts(ps, function(x) x / sum(x) * 100)

# Alpha diversity
alpha_div <- readRDS(file.path(data_dir, "alpha_diversity.rds"))

# Sample metadata
samp_data <- as.data.frame(sample_data(ps))

#-------------------------------------------------------------------------------
# FIGURE 1: ALPHA DIVERSITY PANEL
#-------------------------------------------------------------------------------
cat("\nCreating Figure 1: Alpha Diversity...\n")

create_alpha_panel <- function(data, metrics, comparisons) {
  plots <- lapply(metrics, function(m) {
    ggplot(data, aes(x = group, y = .data[[m$col]], fill = group)) +
      geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.6) +
      geom_jitter(width = 0.15, size = 1.8, alpha = 0.7) +
      scale_fill_manual(values = group_colors) +
      stat_compare_means(
        comparisons = comparisons,
        method = "wilcox.test",
        label = "p.signif",
        size = 3.5,
        step.increase = 0.08
      ) +
      labs(x = "", y = m$label, title = m$title) +
      theme_publication() +
      theme(
        legend.position = "none",
        plot.title = element_text(size = 11, face = "bold"),
        axis.title.y = element_text(size = 10)
      )
  })
  return(plots)
}

comparisons <- list(
  c("Control", "LC_COPD"),
  c("Control", "LC_woCOPD"),
  c("LC_COPD", "LC_woCOPD")
)

metrics <- list(
  list(col = "Observed", label = "Observed ASVs", title = "A) Richness"),
  list(col = "Shannon", label = "Shannon Index", title = "B) Shannon"),
  list(col = "Chao1", label = "Chao1 Estimate", title = "C) Chao1"),
  list(col = "Pielou", label = "Pielou's J", title = "D) Evenness")
)

alpha_panels <- create_alpha_panel(alpha_div, metrics, comparisons)

fig1 <- wrap_plots(alpha_panels, ncol = 2) +
  plot_annotation(
    title = "Figure 1. Alpha Diversity Comparison",
    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0))
  )

save_pub_plot(fig1, "Figure1_alpha_diversity", width = 10, height = 9)

#-------------------------------------------------------------------------------
# FIGURE 2: BETA DIVERSITY
#-------------------------------------------------------------------------------
cat("Creating Figure 2: Beta Diversity...\n")

# Calculate ordinations
ord_bray <- ordinate(ps_rel, method = "PCoA", distance = "bray")
ord_unifrac <- ordinate(ps, method = "PCoA", distance = "wunifrac")

# Get variance explained
var_bray <- round(ord_bray$values$Relative_eig[1:2] * 100, 1)
var_unifrac <- round(ord_unifrac$values$Relative_eig[1:2] * 100, 1)

# Create PCoA plots
pcoa_bray <- plot_ordination(ps_rel, ord_bray, color = "group") +
  geom_point(size = 3, alpha = 0.85) +
  stat_ellipse(aes(color = group),
    type = "norm", level = 0.95,
    linewidth = 0.7, linetype = "dashed"
  ) +
  scale_color_manual(values = group_colors, name = "") +
  labs(
    x = sprintf("PCoA1 (%.1f%%)", var_bray[1]),
    y = sprintf("PCoA2 (%.1f%%)", var_bray[2]),
    title = "A) Bray-Curtis"
  ) +
  theme_publication() +
  theme(legend.position = "bottom")

pcoa_unifrac <- plot_ordination(ps, ord_unifrac, color = "group") +
  geom_point(size = 3, alpha = 0.85) +
  stat_ellipse(aes(color = group),
    type = "norm", level = 0.95,
    linewidth = 0.7, linetype = "dashed"
  ) +
  scale_color_manual(values = group_colors, name = "") +
  labs(
    x = sprintf("PCoA1 (%.1f%%)", var_unifrac[1]),
    y = sprintf("PCoA2 (%.1f%%)", var_unifrac[2]),
    title = "B) Weighted UniFrac"
  ) +
  theme_publication() +
  theme(legend.position = "bottom")

fig2 <- (pcoa_bray | pcoa_unifrac) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Figure 2. Beta Diversity Ordination",
    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0))
  ) &
  theme(legend.position = "bottom")

save_pub_plot(fig2, "Figure2_beta_diversity", width = 12, height = 6)

#-------------------------------------------------------------------------------
# FIGURE 3: TAXONOMIC COMPOSITION
#-------------------------------------------------------------------------------
cat("Creating Figure 3: Taxonomic Composition...\n")

# Phylum-level composition
ps_phylum <- tax_glom(ps_rel, taxrank = "Phylum")
df_phylum <- psmelt(ps_phylum)

# Get top phyla
top_phyla <- df_phylum %>%
  group_by(Phylum) %>%
  summarise(mean_abund = mean(Abundance)) %>%
  arrange(desc(mean_abund)) %>%
  slice_head(n = 8) %>%
  pull(Phylum)

df_phylum <- df_phylum %>%
  mutate(Phylum_Plot = ifelse(Phylum %in% top_phyla, Phylum, "Other")) %>%
  group_by(Sample, group, Phylum_Plot) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")

# Mean by group
phylum_means <- df_phylum %>%
  group_by(group, Phylum_Plot) %>%
  summarise(Mean = mean(Abundance), SE = sd(Abundance) / sqrt(n()), .groups = "drop")

# Order phyla
phylum_order <- phylum_means %>%
  group_by(Phylum_Plot) %>%
  summarise(total = sum(Mean)) %>%
  arrange(desc(total)) %>%
  pull(Phylum_Plot)

phylum_order <- c(setdiff(phylum_order, "Other"), "Other")
phylum_means$Phylum_Plot <- factor(phylum_means$Phylum_Plot, levels = rev(phylum_order))

# Create stacked barplot
phylum_colors <- c(taxa_colors[1:8], "gray70")
names(phylum_colors) <- phylum_order

composition_plot <- ggplot(phylum_means, aes(x = group, y = Mean, fill = Phylum_Plot)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = phylum_colors, name = "Phylum") +
  labs(x = "", y = "Mean Relative Abundance (%)") +
  theme_publication() +
  theme(legend.position = "right") +
  guides(fill = guide_legend(reverse = TRUE))

save_pub_plot(composition_plot, "Figure3_composition", width = 8, height = 6)

#-------------------------------------------------------------------------------
# FIGURE 4: DIFFERENTIAL ABUNDANCE (if ANCOM-BC2 results exist)
#-------------------------------------------------------------------------------
ancom_file <- file.path(data_dir, "ancombc2_output.rds")

if (file.exists(ancom_file)) {
  cat("Creating Figure 4: Differential Abundance...\n")

  ancom_output <- readRDS(ancom_file)

  # Get taxonomy
  tax_df <- as.data.frame(tax_table(tax_glom(ps, taxrank = "Genus", NArm = FALSE)))
  tax_df$taxon <- rownames(tax_df)
  tax_df$Taxon_Label <- ifelse(is.na(tax_df$Genus),
    paste0("Unclassified_", tax_df$Family),
    tax_df$Genus
  )

  res <- ancom_output$res %>%
    left_join(tax_df %>% select(taxon, Taxon_Label), by = "taxon")

  # Filter significant and get top
  sig_data <- res %>%
    filter(diff_groupLC_COPD | diff_groupLC_woCOPD) %>%
    pivot_longer(
      cols = c(lfc_groupLC_COPD, lfc_groupLC_woCOPD),
      names_to = "Comparison",
      values_to = "LFC"
    ) %>%
    mutate(
      Comparison = gsub("lfc_group", "", Comparison),
      Comparison = gsub("_", " ", Comparison)
    ) %>%
    group_by(Taxon_Label) %>%
    filter(abs(LFC) == max(abs(LFC))) %>%
    ungroup() %>%
    arrange(desc(abs(LFC))) %>%
    slice_head(n = 20)

  if (nrow(sig_data) > 0) {
    sig_data$Direction <- ifelse(sig_data$LFC > 0, "Enriched", "Depleted")
    sig_data$Taxon_Label <- factor(sig_data$Taxon_Label,
      levels = sig_data$Taxon_Label[order(sig_data$LFC)]
    )

    da_plot <- ggplot(sig_data, aes(x = LFC, y = Taxon_Label, fill = Direction)) +
      geom_col(width = 0.7, alpha = 0.85) +
      geom_vline(xintercept = 0, color = "gray40", linewidth = 0.5) +
      scale_fill_manual(values = c("Enriched" = "#E64B35", "Depleted" = "#4DBBD5")) +
      facet_wrap(~Comparison, scales = "free_x") +
      labs(
        x = "Log2 Fold Change", y = "",
        title = "Figure 4. Differentially Abundant Taxa (ANCOM-BC2)"
      ) +
      theme_publication() +
      theme(
        legend.position = "top",
        strip.text = element_text(size = 11, face = "bold")
      )

    save_pub_plot(da_plot, "Figure4_differential_abundance", width = 12, height = 8)
  }
}

#-------------------------------------------------------------------------------
# FIGURE 5: FUNCTIONAL PATHWAY ANALYSIS (PCA + Significant Barplot)
#-------------------------------------------------------------------------------
cat("Creating Figure 5: Functional Pathway Analysis...\n")

# Load metadata
metadata <- read.table(metadata_file,
  header = TRUE, sep = "\t",
  row.names = 1, stringsAsFactors = FALSE
)
metadata$group <- factor(metadata$group, levels = c("Control", "LC_COPD", "LC_woCOPD"))

# Load pathway data
pathway_file <- file.path(tables_dir, "pathways_with_descriptions.tsv")
if (!file.exists(pathway_file)) {
  pathway_file <- file.path(project_dir, "results/picrust2/pathways_with_descriptions.tsv")
}

fig5_ok <- FALSE

if (file.exists(pathway_file)) {
  pathways <- read.table(pathway_file,
    header = TRUE, sep = "\t",
    row.names = 1, check.names = FALSE,
    quote = "", comment.char = ""
  )

  # Separate description column if present
  if ("description" %in% colnames(pathways)) {
    pathway_desc <- pathways$description
    names(pathway_desc) <- rownames(pathways)
    pathways <- pathways[, !colnames(pathways) %in% "description"]
  }

  # Match samples
  common_samples <- intersect(colnames(pathways), rownames(metadata))
  if (length(common_samples) > 0) {
    pathways <- pathways[, common_samples]
    meta_func <- metadata[common_samples, , drop = FALSE]

    # Relative abundance
    pathway_rel <- apply(pathways, 2, function(x) x / sum(x) * 100)

    # --- PCA ---
    pathway_t <- t(pathway_rel)
    pca_result <- prcomp(pathway_t, scale. = TRUE, center = TRUE)
    var_explained <- round(pca_result$sdev^2 / sum(pca_result$sdev^2) * 100, 1)

    pca_df <- data.frame(
      PC1 = pca_result$x[, 1],
      PC2 = pca_result$x[, 2],
      group = meta_func$group
    )

    # PERMANOVA
    pathway_dist <- vegdist(pathway_t, method = "bray")
    set.seed(42)
    permanova <- adonis2(pathway_dist ~ group, data = meta_func, permutations = 999)

    pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = group)) +
      geom_point(size = 3.5, alpha = 0.85) +
      stat_ellipse(type = "norm", level = 0.95, linewidth = 0.8, linetype = "dashed") +
      scale_color_manual(values = group_colors) +
      labs(
        x = sprintf("PC1 (%.1f%%)", var_explained[1]),
        y = sprintf("PC2 (%.1f%%)", var_explained[2]),
        title = "A) PCA of Predicted Pathways",
        color = "Group"
      ) +
      annotate("text",
        x = Inf, y = Inf,
        label = sprintf(
          "PERMANOVA\nR\u00b2 = %.3f\np = %.3f",
          permanova$R2[1], permanova$`Pr(>F)`[1]
        ),
        hjust = 1.1, vjust = 1.1, size = 3.5, fontface = "italic"
      ) +
      theme_publication() +
      theme(legend.position = "bottom")

    # --- Kruskal-Wallis testing ---
    p_values <- apply(pathway_rel, 1, function(x) {
      kruskal.test(as.numeric(x) ~ meta_func$group)$p.value
    })

    group_idx <- split(seq_along(meta_func$group), meta_func$group)
    means_control <- rowMeans(pathway_rel[, group_idx[["Control"]], drop = FALSE])
    means_lc_copd <- rowMeans(pathway_rel[, group_idx[["LC_COPD"]], drop = FALSE])
    means_lc_wocopd <- rowMeans(pathway_rel[, group_idx[["LC_woCOPD"]], drop = FALSE])

    kw_results <- data.frame(
      Pathway = rownames(pathway_rel),
      Description = if (exists("pathway_desc")) pathway_desc[rownames(pathway_rel)] else rownames(pathway_rel),
      Mean_Control = means_control,
      Mean_LC_COPD = means_lc_copd,
      Mean_LC_woCOPD = means_lc_wocopd,
      p_value = p_values,
      stringsAsFactors = FALSE
    )
    kw_results$q_value <- p.adjust(kw_results$p_value, method = "BH")
    kw_results <- kw_results %>% arrange(p_value)

    sig_pathways <- kw_results %>% filter(p_value < 0.05)

    # --- Barplot of top significant pathways ---
    if (nrow(sig_pathways) > 0) {
      top_sig <- sig_pathways %>% slice_head(n = 15)
      plot_data <- top_sig %>%
        mutate(Pathway_Short = substr(Description, 1, 45)) %>%
        mutate(Pathway_Short = factor(Pathway_Short, levels = rev(Pathway_Short))) %>%
        pivot_longer(
          cols = starts_with("Mean_"),
          names_to = "Group", values_to = "Abundance"
        ) %>%
        mutate(
          Group = gsub("Mean_", "", Group),
          Group = factor(Group, levels = c("Control", "LC_COPD", "LC_woCOPD"))
        )

      barplot_sig <- ggplot(plot_data, aes(
        x = Pathway_Short,
        y = Abundance, fill = Group
      )) +
        geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
        coord_flip() +
        scale_fill_manual(values = group_colors) +
        labs(
          x = "", y = "Mean Relative Abundance (%)",
          title = "B) Top Differentially Abundant Pathways (p < 0.05)"
        ) +
        theme_publication() +
        theme(
          axis.text.y = element_text(size = 7),
          legend.position = "bottom"
        )

      fig5 <- (pca_plot | barplot_sig) +
        plot_layout(widths = c(1, 1.3), guides = "collect") +
        plot_annotation(
          title = "Figure 5. Functional Pathway Analysis (PICRUSt2)",
          theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0))
        ) &
        theme(legend.position = "bottom")

      fig5_ok <- TRUE
    } else {
      # No significant pathways — show PCA only
      fig5 <- pca_plot +
        plot_annotation(
          title = "Figure 5. Functional Pathway Analysis (PICRUSt2)",
          theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0))
        )
      fig5_ok <- TRUE
    }

    if (fig5_ok) {
      save_pub_plot(fig5, "Figure5_functional_pathways", width = 16, height = 8)
    }
  }
} else {
  cat("  Pathway file not found — skipping Figure 5.\n")
}

#-------------------------------------------------------------------------------
# SUPPLEMENTARY FIGURE S1: RAREFACTION CURVES
#-------------------------------------------------------------------------------
cat("Creating Supplementary Figure S1: Rarefaction...\n")

rarefaction_src <- file.path(output_dir, "alpha_diversity/rarefaction_curves.pdf")
if (file.exists(rarefaction_src)) {
  file.copy(rarefaction_src, file.path(pub_dir, "FigureS1_rarefaction.pdf"),
    overwrite = TRUE
  )
  cat("  Copied rarefaction curves.\n")
}

# Also copy PNG if exists
rarefaction_png <- file.path(output_dir, "alpha_diversity/rarefaction_curves.png")
if (file.exists(rarefaction_png)) {
  file.copy(rarefaction_png, file.path(pub_dir, "FigureS1_rarefaction.png"),
    overwrite = TRUE
  )
}

#-------------------------------------------------------------------------------
# FIGURE S2: VOLCANO PLOTS (COMBINED ANCOM-BC2)
#-------------------------------------------------------------------------------
cat("Creating Figure S2: Volcano Plots...\n")

if (file.exists(ancom_file)) {
  # Re-use ancom_output loaded earlier
  if (!exists("ancom_output")) {
    ancom_output <- readRDS(ancom_file)
  }

  # Build taxonomy reference
  ps_genus_pub <- tax_glom(ps, taxrank = "Genus", NArm = FALSE)
  tax_df_vol <- as.data.frame(tax_table(ps_genus_pub))
  tax_df_vol$taxon <- rownames(tax_df_vol)
  tax_df_vol$Taxon_Label <- ifelse(is.na(tax_df_vol$Genus),
    paste0("Unclassified_", tax_df_vol$Family),
    tax_df_vol$Genus
  )

  res_primary <- ancom_output$res %>%
    left_join(tax_df_vol %>% select(taxon, Taxon_Label, Phylum, Family, Genus),
      by = "taxon"
    )

  # Add p-value significance flags
  add_pval_diff_columns <- function(df) {
    p_cols <- grep("^p_group", colnames(df), value = TRUE)
    for (pc in p_cols) {
      new_col <- sub("^p_", "diff_pval_", pc)
      df[[new_col]] <- !is.na(df[[pc]]) & df[[pc]] < 0.05
    }
    return(df)
  }
  res_primary <- add_pval_diff_columns(res_primary)

  # Volcano plot function (self-contained)
  create_volcano_pub <- function(data, lfc_col, pval_col, sig_col_name, title,
                                 yaxis_label = "-log10(p-value)") {
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

    # Label top taxa — top 8 enriched + top 8 depleted
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
      geom_text_repel(
        data = top_taxa,
        aes(label = Taxon_Label),
        size = 2.8,
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

  # Create three volcano plots
  volcano_lc_copd <- create_volcano_pub(
    res_primary,
    "lfc_groupLC_COPD", "p_groupLC_COPD", "diff_pval_groupLC_COPD",
    "A) LC_COPD vs Control"
  )

  volcano_lc_wocopd <- create_volcano_pub(
    res_primary,
    "lfc_groupLC_woCOPD", "p_groupLC_woCOPD", "diff_pval_groupLC_woCOPD",
    "B) LC_woCOPD vs Control"
  )

  # Try pairwise comparison
  res_pairwise <- ancom_output$res_pair %>%
    left_join(tax_df_vol %>% select(taxon, Taxon_Label, Phylum, Family, Genus),
      by = "taxon"
    )
  res_pairwise <- add_pval_diff_columns(res_pairwise)

  pair_lfc_col <- grep("^lfc_.*LC_COPD.*LC_woCOPD", colnames(res_pairwise), value = TRUE)[1]
  if (is.na(pair_lfc_col)) {
    pair_lfc_col <- grep("^lfc_.*LC_woCOPD.*LC_COPD", colnames(res_pairwise), value = TRUE)[1]
  }
  has_pairwise <- !is.na(pair_lfc_col)

  if (has_pairwise) {
    pair_suffix <- sub("^lfc_", "", pair_lfc_col)
    pair_p_col <- paste0("p_", pair_suffix)
    pair_diff_pval_col <- paste0("diff_pval_", pair_suffix)

    volcano_pair <- create_volcano_pub(
      res_pairwise,
      pair_lfc_col,
      pair_p_col,
      pair_diff_pval_col,
      "C) LC_COPD vs LC_woCOPD"
    )

    fig_s2 <- (volcano_lc_copd | volcano_lc_wocopd | volcano_pair) +
      plot_layout(guides = "collect") +
      plot_annotation(
        title = "Figure S2. Volcano Plots — ANCOM-BC2 Differential Abundance",
        theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0))
      )
    save_pub_plot(fig_s2, "FigureS2_volcano_plots", width = 20, height = 7)
  } else {
    fig_s2 <- (volcano_lc_copd | volcano_lc_wocopd) +
      plot_layout(guides = "collect") +
      plot_annotation(
        title = "Figure S2. Volcano Plots — ANCOM-BC2 Differential Abundance",
        theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0))
      )
    save_pub_plot(fig_s2, "FigureS2_volcano_plots", width = 14, height = 7)
  }
} else {
  cat("  ANCOM-BC2 output not found — skipping Figure S2.\n")
}

#-------------------------------------------------------------------------------
# FIGURE S3: GENUS-LEVEL HEATMAP + CORE MICROBIOME VENN
#-------------------------------------------------------------------------------
cat("Creating Figure S3: Genus Heatmap + Core Venn...\n")

# --- Panel A: Genus-level heatmap (top 30) ---
ps_genus_full_raw <- tax_glom(ps, taxrank = "Genus", NArm = FALSE)
ps_genus_full <- transform_sample_counts(ps_genus_full_raw, function(x) x / sum(x) * 100)
top30_idx <- order(taxa_sums(ps_genus_full), decreasing = TRUE)[1:30]
ps_top30 <- prune_taxa(taxa_names(ps_genus_full)[top30_idx], ps_genus_full)

abund_mat <- as(otu_table(ps_top30), "matrix")
if (!taxa_are_rows(ps_top30)) abund_mat <- t(abund_mat)

# Resolve genus names
tax_top30 <- as.data.frame(tax_table(ps_top30))
raw_genus_names <- ifelse(
  is.na(tax_top30$Genus),
  ifelse(!is.na(tax_top30$Family),
    paste0("Uncl. ", tax_top30$Family),
    ifelse(!is.na(tax_top30$Order),
      paste0("Uncl. ", tax_top30$Order),
      "Unclassified"
    )
  ),
  tax_top30$Genus
)
rownames(abund_mat) <- make.unique(raw_genus_names, sep = "_")

annotation_col <- data.frame(
  Group = samp_data$group,
  row.names = rownames(samp_data)
)
ann_colors <- list(Group = group_colors)

# Save heatmap as PNG at 600 DPI
tryCatch(
  {
    png(file.path(pub_dir, "FigureS3a_genus_heatmap.png"),
      width = 12, height = 10, units = "in", res = PUB_DPI
    )
    pheatmap(
      log10(abund_mat + 0.01),
      annotation_col = annotation_col,
      annotation_colors = ann_colors,
      color = viridis(100),
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      show_colnames = TRUE,
      fontsize_row = 8,
      fontsize_col = 7,
      main = "Figure S3a. Top 30 Genera (log10 Relative Abundance)"
    )
    dev.off()
    cat("  Heatmap saved (FigureS3a_genus_heatmap.png).\n")
  },
  error = function(e) {
    cat(sprintf("  WARNING: Heatmap failed (%s). Skipping.\n", e$message))
    try(dev.off(), silent = TRUE)
  }
)

# Also save as PDF
tryCatch(
  {
    pdf(file.path(pub_dir, "FigureS3a_genus_heatmap.pdf"), width = 12, height = 10)
    pheatmap(
      log10(abund_mat + 0.01),
      annotation_col = annotation_col,
      annotation_colors = ann_colors,
      color = viridis(100),
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      show_colnames = TRUE,
      fontsize_row = 8,
      fontsize_col = 7,
      main = "Figure S3a. Top 30 Genera (log10 Relative Abundance)"
    )
    dev.off()
  },
  error = function(e) {
    try(dev.off(), silent = TRUE)
  }
)

# --- Panel B: Core Microbiome Venn ---
# Copy the existing Venn diagram
venn_src <- file.path(output_dir, "composition/core_microbiome_venn.png")
if (file.exists(venn_src)) {
  file.copy(venn_src, file.path(pub_dir, "FigureS3b_core_venn.png"), overwrite = TRUE)
  cat("  Copied Venn diagram (FigureS3b_core_venn.png).\n")
}

# If VennDiagram is available, also regenerate at higher resolution
if (requireNamespace("VennDiagram", quietly = TRUE) &&
  requireNamespace("microbiome", quietly = TRUE)) {
  library(VennDiagram)
  library(microbiome)

  core_by_group <- lapply(levels(samp_data$group), function(g) {
    ps_group <- prune_samples(sample_data(ps_rel)$group == g, ps_rel)
    core_members(ps_group, detection = 0.1, prevalence = 0.5)
  })
  names(core_by_group) <- levels(samp_data$group)

  venn_data <- list(
    Control = core_by_group$Control,
    LC_COPD = core_by_group$LC_COPD,
    LC_woCOPD = core_by_group$LC_woCOPD
  )

  tryCatch(
    {
      venn.diagram(
        x = venn_data,
        category.names = names(venn_data),
        filename = file.path(pub_dir, "FigureS3b_core_venn_hires.png"),
        output = TRUE,
        imagetype = "png",
        height = 3000,
        width = 3000,
        resolution = PUB_DPI,
        fill = unname(group_colors[names(venn_data)]),
        alpha = 0.5,
        cat.cex = 1.4,
        cex = 1.8,
        main = "Core Microbiome Overlap",
        main.cex = 1.6
      )
      cat("  High-resolution Venn diagram saved.\n")
    },
    error = function(e) {
      cat(sprintf("  WARNING: Venn diagram failed (%s).\n", e$message))
    }
  )
}

#-------------------------------------------------------------------------------
# FIGURE S4: KO & EC SIGNIFICANT BARPLOTS
#-------------------------------------------------------------------------------
cat("Creating Figure S4: KO & EC Barplots...\n")

# Helper: load functional data, run KW, and create barplot
create_functional_barplot <- function(data_file, meta, id_col, title_label) {
  if (!file.exists(data_file)) {
    cat(sprintf("  %s file not found — skipping.\n", id_col))
    return(NULL)
  }

  func_data <- read.table(data_file,
    header = TRUE, sep = "\t",
    row.names = 1, check.names = FALSE,
    quote = "", comment.char = ""
  )

  # Separate description
  func_desc <- NULL
  if ("description" %in% colnames(func_data)) {
    func_desc <- func_data$description
    names(func_desc) <- rownames(func_data)
    func_data <- func_data[, !colnames(func_data) %in% "description"]
  }

  # Match samples
  common <- intersect(colnames(func_data), rownames(meta))
  if (length(common) == 0) {
    return(NULL)
  }

  func_data <- func_data[, common]
  meta_sub <- meta[common, , drop = FALSE]

  # Relative abundance
  func_rel <- apply(func_data, 2, function(x) x / sum(x) * 100)

  # Remove zero-variance rows
  row_var <- apply(func_rel, 1, var)
  func_rel <- func_rel[row_var > 0, ]

  # Kruskal-Wallis
  p_values <- apply(func_rel, 1, function(x) {
    kruskal.test(as.numeric(x) ~ meta_sub$group)$p.value
  })

  group_idx <- split(seq_along(meta_sub$group), meta_sub$group)
  means_ctrl <- rowMeans(func_rel[, group_idx[["Control"]], drop = FALSE])
  means_copd <- rowMeans(func_rel[, group_idx[["LC_COPD"]], drop = FALSE])
  means_wocopd <- rowMeans(func_rel[, group_idx[["LC_woCOPD"]], drop = FALSE])

  results <- data.frame(
    ID = rownames(func_rel),
    Description = if (!is.null(func_desc)) func_desc[rownames(func_rel)] else rownames(func_rel),
    Mean_Control = means_ctrl,
    Mean_LC_COPD = means_copd,
    Mean_LC_woCOPD = means_wocopd,
    p_value = p_values,
    stringsAsFactors = FALSE
  )
  results$q_value <- p.adjust(results$p_value, method = "BH")
  results <- results %>% arrange(p_value)

  sig_results <- results %>% filter(p_value < 0.05)

  if (nrow(sig_results) == 0) {
    cat(sprintf("  No significant %s features.\n", id_col))
    return(NULL)
  }

  top_sig <- sig_results %>% slice_head(n = 15)

  # Create labels: ID + description (truncated)
  top_sig$Label <- ifelse(
    !is.na(top_sig$Description) & top_sig$Description != top_sig$ID,
    paste0(top_sig$ID, ": ", substr(top_sig$Description, 1, 35)),
    top_sig$ID
  )

  plot_data <- top_sig %>%
    mutate(Label = factor(Label, levels = rev(Label))) %>%
    pivot_longer(
      cols = starts_with("Mean_"),
      names_to = "Group", values_to = "Abundance"
    ) %>%
    mutate(
      Group = gsub("Mean_", "", Group),
      Group = factor(Group, levels = c("Control", "LC_COPD", "LC_woCOPD"))
    )

  p <- ggplot(plot_data, aes(x = Label, y = Abundance, fill = Group)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    coord_flip() +
    scale_fill_manual(values = group_colors) +
    labs(
      x = "", y = "Mean Relative Abundance (%)",
      title = title_label
    ) +
    theme_publication() +
    theme(
      axis.text.y = element_text(size = 6.5),
      legend.position = "bottom"
    )

  return(p)
}

# Load metadata for functional analysis
meta_func <- read.table(metadata_file,
  header = TRUE, sep = "\t",
  row.names = 1, stringsAsFactors = FALSE
)
meta_func$group <- factor(meta_func$group, levels = c("Control", "LC_COPD", "LC_woCOPD"))

# KO barplot
ko_file <- file.path(tables_dir, "KO_metagenome_with_descriptions.tsv")
ko_plot <- create_functional_barplot(
  ko_file, meta_func, "KO",
  "A) Top Differentially Abundant KOs (p < 0.05)"
)

# EC barplot
ec_file <- file.path(tables_dir, "EC_metagenome_with_descriptions.tsv")
ec_plot <- create_functional_barplot(
  ec_file, meta_func, "EC",
  "B) Top Differentially Abundant ECs (p < 0.05)"
)

# Combine KO + EC
if (!is.null(ko_plot) && !is.null(ec_plot)) {
  fig_s4 <- (ko_plot / ec_plot) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = "Figure S4. Functional Analysis — KO & EC Significant Features",
      theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0))
    ) &
    theme(legend.position = "bottom")

  save_pub_plot(fig_s4, "FigureS4_ko_ec_barplots", width = 14, height = 16)
} else if (!is.null(ko_plot)) {
  save_pub_plot(ko_plot, "FigureS4_ko_barplot", width = 14, height = 8)
} else if (!is.null(ec_plot)) {
  save_pub_plot(ec_plot, "FigureS4_ec_barplot", width = 14, height = 8)
} else {
  cat("  Neither KO nor EC data available — skipping Figure S4.\n")
}

#-------------------------------------------------------------------------------
# SUMMARY STATISTICS TABLE
#-------------------------------------------------------------------------------
cat("\nGenerating summary statistics table...\n")

# Sample size and sequencing depth
summary_stats <- samp_data %>%
  group_by(group) %>%
  summarise(
    n = n(),
    .groups = "drop"
  )

# Add read counts
read_counts <- read.csv(file.path(data_dir, "sample_read_counts.csv"))
read_summary <- read_counts %>%
  group_by(group) %>%
  summarise(
    Mean_Reads = round(mean(reads)),
    SD_Reads = round(sd(reads)),
    Min_Reads = min(reads),
    Max_Reads = max(reads),
    .groups = "drop"
  )

summary_stats <- left_join(summary_stats, read_summary, by = "group")

# Add ASV counts
asv_summary <- alpha_div %>%
  group_by(group) %>%
  summarise(
    Mean_ASVs = round(mean(Observed)),
    SD_ASVs = round(sd(Observed)),
    .groups = "drop"
  )

summary_stats <- left_join(summary_stats, asv_summary, by = "group")

write.csv(summary_stats, file.path(pub_dir, "Table1_sample_summary.csv"),
  row.names = FALSE
)

cat("\nSample Summary:\n")
print(summary_stats)

#-------------------------------------------------------------------------------
# CLEAN UP AND FINISH
#-------------------------------------------------------------------------------
cat("\n============================================\n")
cat("Publication Figures Complete!\n")
cat("============================================\n")
cat("\nGenerated files in:", pub_dir, "\n")
cat("\nAll PNG figures saved at", PUB_DPI, "DPI\n\n")

pub_files <- list.files(pub_dir, pattern = "\\.(pdf|png|tiff|csv)$")
cat(sprintf("Total files: %d\n", length(pub_files)))
for (f in pub_files) {
  cat(sprintf("  %s\n", f))
}
