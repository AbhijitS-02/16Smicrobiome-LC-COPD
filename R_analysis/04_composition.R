#===============================================================================
# 04_COMPOSITION.R - Taxonomic Composition Analysis
#===============================================================================
# Relative abundance plots, heatmaps, and core microbiome analysis
#===============================================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(microbiome)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(viridis)
  library(pheatmap)
  library(RColorBrewer)
})

source("D:/LC_COPD_microbiome/R_analysis/functions/plotting_theme.R")

# Configuration
project_dir <- "D:/LC_COPD_microbiome"
data_dir <- file.path(project_dir, "R_analysis/data")
output_dir <- file.path(project_dir, "figures/composition")
tables_dir <- file.path(project_dir, "results/tables")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("============================================\n")
cat("Taxonomic Composition Analysis\n")
cat("============================================\n\n")

#-------------------------------------------------------------------------------
# LOAD DATA
#-------------------------------------------------------------------------------
ps <- readRDS(file.path(data_dir, "phyloseq_silva_filtered.rds"))
ps_rel <- transform_sample_counts(ps, function(x) x / sum(x) * 100)

samp_data <- as.data.frame(sample_data(ps))

#-------------------------------------------------------------------------------
# PHYLUM-LEVEL BARPLOT
#-------------------------------------------------------------------------------
cat("Creating phylum-level barplot...\n")

# Aggregate to phylum on RAW COUNTS first, then convert to relative abundance
# (glomming integers is faster and avoids floating-point summation drift)
ps_phylum_raw <- tax_glom(ps, taxrank = "Phylum")
ps_phylum <- transform_sample_counts(ps_phylum_raw, function(x) x / sum(x) * 100)

# Get top phyla
phylum_sums <- taxa_sums(ps_phylum)
top_phyla <- names(sort(phylum_sums, decreasing = TRUE))[1:10]

# Melt data
df_phylum <- psmelt(ps_phylum) %>%
  mutate(Phylum = ifelse(OTU %in% top_phyla, Phylum, "Other")) %>%
  group_by(Sample, group, Phylum) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")

# Order phyla by mean abundance
phylum_order <- df_phylum %>%
  group_by(Phylum) %>%
  summarise(mean_abund = mean(Abundance)) %>%
  arrange(desc(mean_abund)) %>%
  pull(Phylum)

# Move "Other" to end
phylum_order <- c(setdiff(phylum_order, "Other"), "Other")
df_phylum$Phylum <- factor(df_phylum$Phylum, levels = rev(phylum_order))

# Vibrant phylum palette — high-contrast, matching genus palette style
phylum_base_palette <- c(
  "#845EF7", "#51CF66", "#FF922B", "#339AF0", "#FF6B6B",
  "#20C997", "#CC5DE8", "#FCC419", "#F06595", "#4DABF7",
  "#94D82D", "#DA77F2", "#FD7E14", "#66D9E8", "#FAB005"
)
n_phyla <- length(unique(df_phylum$Phylum))
# "Other" is the first level after rev(); assign it a neutral gray
phylum_colors <- c("#E8E8E8", phylum_base_palette[1:(n_phyla - 1)])
names(phylum_colors) <- levels(df_phylum$Phylum)

# Order samples within each group by ascending "Other" abundance
phylum_sample_order <- df_phylum %>%
  filter(Phylum == "Other") %>%
  arrange(group, Abundance) %>%
  pull(Sample)
# Samples with zero "Other" won't appear above; append them at the front
all_phylum_samples <- unique(df_phylum$Sample)
phylum_sample_order <- c(setdiff(all_phylum_samples, phylum_sample_order), phylum_sample_order)
df_phylum$Sample <- factor(df_phylum$Sample, levels = phylum_sample_order)

# Plot
phylum_plot <- ggplot(df_phylum, aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", width = 0.95) +
  facet_grid(~ group, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = phylum_colors) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0))) +
  labs(x = "", y = "Relative Abundance (%)", 
       title = "Phylum-Level Composition") +
  theme_publication() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    legend.position = "right",
    strip.text = element_text(face = "bold", size = 12, margin = margin(t = 6, b = 6)),
    strip.background = element_rect(fill = "gray95", color = "black", linewidth = 0.5)
  ) +
  guides(fill = guide_legend(reverse = TRUE, ncol = 1))

ggsave(file.path(output_dir, "barplot_phylum.pdf"), 
       phylum_plot, width = 14, height = 6, dpi = 300)
ggsave(file.path(output_dir, "barplot_phylum.png"), 
       phylum_plot, width = 14, height = 6, dpi = 300)

#-------------------------------------------------------------------------------
# GENUS-LEVEL BARPLOT (TOP 30)
#-------------------------------------------------------------------------------
cat("Creating genus-level barplot...\n")

# Aggregate to genus on RAW COUNTS first, then convert to relative abundance
ps_genus_raw <- tax_glom(ps, taxrank = "Genus", NArm = FALSE)
ps_genus <- transform_sample_counts(ps_genus_raw, function(x) x / sum(x) * 100)

# Get top genera
genus_sums <- taxa_sums(ps_genus)
top_genera_idx <- order(genus_sums, decreasing = TRUE)[1:30]
top_genera <- taxa_names(ps_genus)[top_genera_idx]

# Resolve unclassified genera using the lowest available taxonomic rank
# (e.g. "Uncl. Lachnospiraceae" instead of a generic "Unclassified")
genus_tax_df <- as.data.frame(tax_table(ps_genus))
resolved_genus <- ifelse(
  is.na(genus_tax_df$Genus),
  ifelse(!is.na(genus_tax_df$Family),
         paste0("Uncl. ", genus_tax_df$Family),
         ifelse(!is.na(genus_tax_df$Order),
                paste0("Uncl. ", genus_tax_df$Order),
                ifelse(!is.na(genus_tax_df$Class),
                       paste0("Uncl. ", genus_tax_df$Class),
                       ifelse(!is.na(genus_tax_df$Phylum),
                              paste0("Uncl. ", genus_tax_df$Phylum),
                              "Unclassified")))),
  genus_tax_df$Genus
)
# Ensure uniqueness across ASVs that share the same resolved name
resolved_genus <- make.unique(resolved_genus, sep = " ")
# Map resolved names back to ASV IDs for downstream labelling
names(resolved_genus) <- taxa_names(ps_genus)

# Melt and process
df_genus <- psmelt(ps_genus) %>%
  mutate(Genus = resolved_genus[OTU]) %>%
  mutate(Genus = ifelse(OTU %in% top_genera, Genus, "Other")) %>%
  group_by(Sample, group, Genus) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")

# Order by abundance
genus_order <- df_genus %>%
  group_by(Genus) %>%
  summarise(mean_abund = mean(Abundance)) %>%
  arrange(desc(mean_abund)) %>%
  pull(Genus)

# Only append special categories that actually exist in the data
existing_genera <- unique(df_genus$Genus)
special_tail <- c()
if ("Unclassified" %in% existing_genera) special_tail <- c(special_tail, "Unclassified")
if ("Other" %in% existing_genera) special_tail <- c(special_tail, "Other")

genus_order <- c(setdiff(genus_order, c("Other", "Unclassified")), special_tail)
df_genus$Genus <- factor(df_genus$Genus, levels = rev(genus_order))

# Vibrant genus palette — 30 high-contrast colors (reversed)
genus_base_hues <- c(
  "#845EF7", "#94D82D", "#20C997", "#FD7E14", "#CC5DE8",
  "#339AF0", "#51CF66", "#FCC419", "#5C7CFA", "#FF922B",
  "#C0EB75", "#91A7FF", "#FFA8A8", "#748FFC", "#38D9A9",
  "#E599F7", "#FAB005", "#A9E34B", "#74C0FC", "#FF8787",
  "#DA77F2", "#63E6BE", "#66D9E8", "#F06595", "#B197FC",
  "#4DABF7", "#69DB7C", "#FFD43B", "#FFA94D", "#FF6B6B"
)
# Derive count from the finalized factor levels so colors always match
n_special <- length(special_tail)
n_genera <- length(levels(df_genus$Genus))
n_regular <- n_genera - n_special
# Grays go FIRST — after rev(), special categories ("Other" etc.) are the first levels
special_grays <- c("#E8E8E8", "#D5D5D5")[seq_len(n_special)]
genus_pal <- if (n_regular <= length(genus_base_hues)) {
  genus_base_hues[1:n_regular]
} else {
  colorRampPalette(genus_base_hues)(n_regular)
}
genus_colors <- c(special_grays, genus_pal)
names(genus_colors) <- levels(df_genus$Genus)

# Order samples within each group by ascending "Other" abundance
genus_sample_order <- df_genus %>%
  filter(Genus == "Other") %>%
  group_by(Sample, group) %>%
  summarise(unclass_total = sum(Abundance), .groups = "drop") %>%
  arrange(group, unclass_total) %>%
  pull(Sample)
all_genus_samples <- unique(df_genus$Sample)
genus_sample_order <- c(setdiff(all_genus_samples, genus_sample_order), genus_sample_order)
df_genus$Sample <- factor(df_genus$Sample, levels = genus_sample_order)

# Plot
genus_plot <- ggplot(df_genus, aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", width = 0.95) +
  facet_grid(~ group, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = genus_colors) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0))) +
  labs(x = "", y = "Relative Abundance (%)", 
       title = "Genus-Level Composition (Top 30)") +
  theme_publication() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    legend.position = "right",
    strip.text = element_text(face = "bold", size = 12, margin = margin(t = 6, b = 6)),
    strip.background = element_rect(fill = "gray95", color = "black", linewidth = 0.5)
  ) +
  guides(fill = guide_legend(reverse = TRUE, ncol = 1))

ggsave(file.path(output_dir, "barplot_genus.pdf"), 
       genus_plot, width = 14, height = 7, dpi = 300)
ggsave(file.path(output_dir, "barplot_genus.png"), 
       genus_plot, width = 14, height = 7, dpi = 300)

#-------------------------------------------------------------------------------
# MEAN ABUNDANCE BY GROUP (STACKED BAR)
#-------------------------------------------------------------------------------
cat("Creating group mean abundance plots...\n")

# Phylum means by group
phylum_means <- df_phylum %>%
  group_by(group, Phylum) %>%
  summarise(Mean_Abundance = mean(Abundance), .groups = "drop")

group_phylum_plot <- ggplot(phylum_means, aes(x = group, y = Mean_Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = phylum_colors) +
  labs(x = "", y = "Mean Relative Abundance (%)",
       title = "Phylum Composition by Group") +
  theme_publication() +
  theme(legend.position = "right") +
  guides(fill = guide_legend(reverse = TRUE))

ggsave(file.path(output_dir, "barplot_phylum_by_group.pdf"), 
       group_phylum_plot, width = 8, height = 6, dpi = 300)

#-------------------------------------------------------------------------------
# HEATMAP (TOP 30 GENERA)
#-------------------------------------------------------------------------------
cat("Creating genus heatmap...\n")

# Get top 30 genera — glom from raw counts, then transform
ps_genus_full_raw <- tax_glom(ps, taxrank = "Genus", NArm = FALSE)
ps_genus_full <- transform_sample_counts(ps_genus_full_raw, function(x) x / sum(x) * 100)
top30_idx <- order(taxa_sums(ps_genus_full), decreasing = TRUE)[1:30]
ps_top30 <- prune_taxa(taxa_names(ps_genus_full)[top30_idx], ps_genus_full)

# Create abundance matrix
abund_mat <- as(otu_table(ps_top30), "matrix")
if (!taxa_are_rows(ps_top30)) abund_mat <- t(abund_mat)

# Resolve genus names — fall back to the lowest known rank for NAs
tax_top30 <- as.data.frame(tax_table(ps_top30))
raw_genus_names <- ifelse(
  is.na(tax_top30$Genus),
  ifelse(!is.na(tax_top30$Family),
         paste0("Uncl. ", tax_top30$Family),
         ifelse(!is.na(tax_top30$Order),
                paste0("Uncl. ", tax_top30$Order),
                "Unclassified")),
  tax_top30$Genus
)
# Force unique row names (e.g. "Uncl. Lachnospiraceae", "Uncl. Lachnospiraceae_1")
rownames(abund_mat) <- make.unique(raw_genus_names, sep = "_")

# Annotation
annotation_col <- data.frame(
  Group = samp_data$group,
  row.names = rownames(samp_data)
)

ann_colors <- list(Group = group_colors)

# Create heatmap
pdf(file.path(output_dir, "heatmap_top30_genera.pdf"), width = 12, height = 10)
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
  main = "Top 30 Genera (log10 Relative Abundance)"
)
dev.off()

# Save as PNG too
png(file.path(output_dir, "heatmap_top30_genera.png"), 
    width = 12, height = 10, units = "in", res = 300)
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
  main = "Top 30 Genera (log10 Relative Abundance)"
)
dev.off()

#-------------------------------------------------------------------------------
# CORE MICROBIOME ANALYSIS
#-------------------------------------------------------------------------------
cat("Analyzing core microbiome...\n")

# Define core: present in >50% samples with >0.1% abundance
core_taxa <- core_members(ps_rel, detection = 0.1, prevalence = 0.5)

cat(sprintf("  Core taxa (>50%% prevalence, >0.1%% abundance): %d\n", length(core_taxa)))

# Core by group
core_by_group <- lapply(levels(samp_data$group), function(g) {
  ps_group <- prune_samples(sample_data(ps_rel)$group == g, ps_rel)
  core_members(ps_group, detection = 0.1, prevalence = 0.5)
})
names(core_by_group) <- levels(samp_data$group)

cat("  Core taxa per group:\n")
for (g in names(core_by_group)) {
  cat(sprintf("    %s: %d taxa\n", g, length(core_by_group[[g]])))
}

# Create Venn diagram data
venn_data <- list(
  Control = core_by_group$Control,
  LC_COPD = core_by_group$LC_COPD,
  LC_woCOPD = core_by_group$LC_woCOPD
)

# If VennDiagram available, create plot
if (requireNamespace("VennDiagram", quietly = TRUE)) {
  library(VennDiagram)
  
  venn.diagram(
    x = venn_data,
    category.names = names(venn_data),
    filename = file.path(output_dir, "core_microbiome_venn.png"),
    output = TRUE,
    imagetype = "png",
    height = 2000,
    width = 2000,
    resolution = 300,
    fill = unname(group_colors[names(venn_data)]),
    alpha = 0.5,
    cat.cex = 1.2,
    cex = 1.5
  )
}

# Save Venn diagram summary as a text file
cat("  Saving Venn diagram summary...\n")
venn_txt <- file.path(tables_dir, "core_microbiome_venn_summary.txt")
tax_all <- as.data.frame(tax_table(ps))

# Helper: resolve ASV IDs to readable taxonomy strings
resolve_taxonomy <- function(asv_ids, tax_df) {
  sapply(asv_ids, function(id) {
    row <- tax_df[id, ]
    ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
    ranks <- ranks[ranks %in% colnames(row)]
    paste(na.omit(unlist(row[ranks])), collapse = "; ")
  })
}

sink(venn_txt)
cat("================================================================\n")
cat("Core Microbiome Venn Diagram Summary\n")
cat("================================================================\n")
cat(sprintf("Detection threshold: >0.1%% abundance\n"))
cat(sprintf("Prevalence threshold: >50%% of samples\n\n"))

# Per-group core taxa
for (g in names(core_by_group)) {
  taxa_ids <- core_by_group[[g]]
  cat(sprintf("--- %s (%d core taxa) ---\n", g, length(taxa_ids)))
  if (length(taxa_ids) > 0) {
    tax_strings <- resolve_taxonomy(taxa_ids, tax_all)
    for (i in seq_along(taxa_ids)) {
      cat(sprintf("  %s : %s\n", taxa_ids[i], tax_strings[i]))
    }
  }
  cat("\n")
}

# Shared among ALL groups (center of Venn)
shared_all <- Reduce(intersect, core_by_group)
cat(sprintf("--- Shared across ALL groups (%d taxa) ---\n", length(shared_all)))
if (length(shared_all) > 0) {
  tax_strings <- resolve_taxonomy(shared_all, tax_all)
  for (i in seq_along(shared_all)) {
    cat(sprintf("  %s : %s\n", shared_all[i], tax_strings[i]))
  }
} else {
  cat("  (none)\n")
}
cat("\n")

# Pairwise shared
group_names <- names(core_by_group)
pairs <- combn(group_names, 2, simplify = FALSE)
for (pair in pairs) {
  shared <- intersect(core_by_group[[pair[1]]], core_by_group[[pair[2]]])
  cat(sprintf("--- Shared: %s & %s (%d taxa) ---\n", pair[1], pair[2], length(shared)))
  if (length(shared) > 0) {
    tax_strings <- resolve_taxonomy(shared, tax_all)
    for (i in seq_along(shared)) {
      cat(sprintf("  %s : %s\n", shared[i], tax_strings[i]))
    }
  } else {
    cat("  (none)\n")
  }
  cat("\n")
}

# Unique to each group
for (g in group_names) {
  others <- setdiff(group_names, g)
  unique_taxa <- setdiff(core_by_group[[g]], Reduce(union, core_by_group[others]))
  cat(sprintf("--- Unique to %s (%d taxa) ---\n", g, length(unique_taxa)))
  if (length(unique_taxa) > 0) {
    tax_strings <- resolve_taxonomy(unique_taxa, tax_all)
    for (i in seq_along(unique_taxa)) {
      cat(sprintf("  %s : %s\n", unique_taxa[i], tax_strings[i]))
    }
  } else {
    cat("  (none)\n")
  }
  cat("\n")
}

cat("================================================================\n")
sink()
cat(sprintf("  Saved: %s\n", venn_txt))

#-------------------------------------------------------------------------------
# PREVALENCE-ABUNDANCE PLOT
#-------------------------------------------------------------------------------
cat("Creating prevalence-abundance plot...\n")

# Calculate prevalence and mean abundance
prevalence <- apply(otu_table(ps), 1, function(x) sum(x > 0) / length(x) * 100)
mean_abund <- taxa_sums(ps_rel) / nsamples(ps_rel)

prev_abund_df <- data.frame(
  Prevalence = prevalence,
  Abundance = mean_abund,
  Taxon = taxa_names(ps)
)

# Add taxonomy
tax_df <- as.data.frame(tax_table(ps))
prev_abund_df$Phylum <- tax_df$Phylum

prev_abund_plot <- ggplot(prev_abund_df, aes(x = Prevalence, y = Abundance, color = Phylum)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_y_log10() +
  scale_color_manual(values = colorRampPalette(brewer.pal(12, "Set3"))(
    length(unique(prev_abund_df$Phylum)))) +
  labs(x = "Prevalence (%)", y = "Mean Relative Abundance (%, log scale)",
       title = "Prevalence vs Abundance") +
  theme_publication() +
  theme(legend.position = "right")

ggsave(file.path(output_dir, "prevalence_abundance.pdf"), 
       prev_abund_plot, width = 10, height = 7, dpi = 300)

#-------------------------------------------------------------------------------
# SAVE RESULTS
#-------------------------------------------------------------------------------
cat("\nSaving composition tables...\n")

# Save phylum-level table
phylum_table <- df_phylum %>%
  pivot_wider(names_from = Sample, values_from = Abundance, values_fill = 0)
write.csv(phylum_table, file.path(tables_dir, "composition_phylum.csv"), 
          row.names = FALSE)

# Save genus-level table
genus_table <- df_genus %>%
  pivot_wider(names_from = Sample, values_from = Abundance, values_fill = 0)
write.csv(genus_table, file.path(tables_dir, "composition_genus.csv"), 
          row.names = FALSE)

# Save core microbiome
core_df <- data.frame(
  ASV = core_taxa,
  Taxonomy = apply(tax_table(ps)[core_taxa, ], 1, paste, collapse = "; ")
)
write.csv(core_df, file.path(tables_dir, "core_microbiome.csv"), 
          row.names = FALSE)

cat("\n============================================\n")
cat("Composition Analysis Complete!\n")
cat("============================================\n")
