#===============================================================================
# 02_ALPHA_DIVERSITY.R - Alpha Diversity Analysis
#===============================================================================
# Calculate and visualize alpha diversity metrics with statistical comparisons
#===============================================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(vegan)
  library(ggplot2)
  library(ggpubr)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(microbiome)
})

# Source custom plotting functions
source("D:/LC_COPD_microbiome/R_analysis/functions/plotting_theme.R")

# Configuration
project_dir <- "D:/LC_COPD_microbiome"
data_dir <- file.path(project_dir, "R_analysis/data")
output_dir <- file.path(project_dir, "figures/alpha_diversity")
tables_dir <- file.path(project_dir, "results/tables")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("============================================\n")
cat("Alpha Diversity Analysis\n")
cat("============================================\n\n")

#-------------------------------------------------------------------------------
# LOAD DATA
#-------------------------------------------------------------------------------
cat("Loading phyloseq object...\n")
ps <- readRDS(file.path(data_dir, "phyloseq_silva_filtered.rds"))

cat(sprintf("  Taxa: %d\n", ntaxa(ps)))
cat(sprintf("  Samples: %d\n", nsamples(ps)))

#-------------------------------------------------------------------------------
# CALCULATE ALPHA DIVERSITY METRICS
#-------------------------------------------------------------------------------
cat("\nCalculating alpha diversity metrics...\n")

# Get sample data
samp_data <- as.data.frame(sample_data(ps))

# Calculate multiple metrics
# Observed richness, Chao1, ACE
richness <- estimate_richness(ps, measures = c("Observed", "Chao1", "ACE", 
                                                "Shannon", "Simpson", "InvSimpson", 
                                                "Fisher"))

# Faith's phylogenetic diversity (requires tree)
if (!is.null(phy_tree(ps))) {
  cat("  Calculating Faith's PD...\n")
  faith_pd <- microbiome::alpha(ps, index = "diversity_faith")
  richness$Faith_PD <- faith_pd$diversity_faith
}

# Pielou's evenness
cat("  Calculating Pielou's evenness...\n")
evenness <- microbiome::evenness(ps, index = "pielou")
richness$Pielou <- evenness$pielou

# Combine with metadata
alpha_div <- cbind(samp_data, richness)
alpha_div$sample_id <- rownames(alpha_div)

cat("  Calculated metrics:", paste(names(richness), collapse = ", "), "\n")

#-------------------------------------------------------------------------------
# STATISTICAL TESTS
#-------------------------------------------------------------------------------
cat("\nRunning statistical tests...\n")

metrics <- c("Observed", "Shannon", "Simpson", "InvSimpson", "Chao1", "ACE", 
             "Fisher", "Faith_PD", "Pielou")

# Keep only available metrics
metrics <- metrics[metrics %in% names(alpha_div)]

# Kruskal-Wallis tests
kw_results <- lapply(metrics, function(m) {
  test <- kruskal.test(as.formula(paste(m, "~ group")), data = alpha_div)
  data.frame(
    Metric = m,
    Statistic = test$statistic,
    df = test$parameter,
    p_value = test$p.value,
    stringsAsFactors = FALSE
  )
}) %>% bind_rows()

# Adjust p-values
kw_results$p_adjusted <- p.adjust(kw_results$p_value, method = "BH")
kw_results$Significance <- ifelse(kw_results$p_adjusted < 0.001, "***",
                                   ifelse(kw_results$p_adjusted < 0.01, "**",
                                          ifelse(kw_results$p_adjusted < 0.05, "*", "ns")))

cat("\nKruskal-Wallis test results:\n")
print(kw_results)

# Pairwise Wilcoxon tests for significant metrics
pairwise_results <- list()

for (m in metrics[kw_results$p_adjusted[match(metrics, kw_results$Metric)] < 0.05]) {
  pw <- pairwise.wilcox.test(alpha_div[[m]], alpha_div$group, 
                              p.adjust.method = "BH")
  pairwise_results[[m]] <- as.data.frame(pw$p.value) %>%
    rownames_to_column("Group1") %>%
    pivot_longer(-Group1, names_to = "Group2", values_to = "p_adjusted") %>%
    mutate(Metric = m, .before = 1) %>%
    filter(!is.na(p_adjusted))
}

if (length(pairwise_results) > 0) {
  pairwise_df <- bind_rows(pairwise_results)
  cat("\nSignificant pairwise comparisons:\n")
  print(filter(pairwise_df, p_adjusted < 0.05))
}

#-------------------------------------------------------------------------------
# VISUALIZATION: INDIVIDUAL METRIC PLOTS
#-------------------------------------------------------------------------------
cat("\nGenerating visualizations...\n")

# Define comparisons for significance annotation
comparisons <- list(
  c("Control", "LC_COPD"),
  c("Control", "LC_woCOPD"),
  c("LC_COPD", "LC_woCOPD")
)

# Function to create alpha diversity boxplot
plot_alpha_metric <- function(data, metric, ylab = NULL) {
  if (is.null(ylab)) ylab <- metric
  
  p <- ggplot(data, aes(x = group, y = .data[[metric]], fill = group)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 2, alpha = 0.6) +
    scale_fill_manual(values = group_colors) +
    stat_compare_means(comparisons = comparisons, 
                       method = "wilcox.test",
                       label = "p.signif",
                       size = 4,
                       step.increase = 0.1) +
    stat_compare_means(method = "kruskal.test", 
                       label.y.npc = 0.95,
                       label = "p.format",
                       size = 3.5) +
    labs(x = "", y = ylab, title = metric) +
    theme_publication() +
    theme(legend.position = "none")
  
  return(p)
}

# Create individual plots
alpha_plots <- list()

alpha_plots$observed <- plot_alpha_metric(alpha_div, "Observed", "Observed ASVs")
alpha_plots$shannon <- plot_alpha_metric(alpha_div, "Shannon", "Shannon Index")
alpha_plots$simpson <- plot_alpha_metric(alpha_div, "Simpson", "Simpson Index")
alpha_plots$invsimpson <- plot_alpha_metric(alpha_div, "InvSimpson", "Inverse Simpson")
alpha_plots$chao1 <- plot_alpha_metric(alpha_div, "Chao1", "Chao1 Richness")
alpha_plots$ace <- plot_alpha_metric(alpha_div, "ACE", "ACE Richness")

if ("Faith_PD" %in% names(alpha_div)) {
  alpha_plots$faith <- plot_alpha_metric(alpha_div, "Faith_PD", "Faith's PD")
}
alpha_plots$pielou <- plot_alpha_metric(alpha_div, "Pielou", "Pielou's Evenness")

if ("Fisher" %in% names(alpha_div)) {
  alpha_plots$fisher <- plot_alpha_metric(alpha_div, "Fisher", "Fisher's Alpha")
}

#-------------------------------------------------------------------------------
# COMBINED PANEL FIGURE
#-------------------------------------------------------------------------------
cat("Creating combined panel figure...\n")

# Main figure (4 key metrics)
main_panel <- (alpha_plots$observed + alpha_plots$shannon) /
              (alpha_plots$chao1 + alpha_plots$pielou) +
  plot_annotation(
    title = "Alpha Diversity",
    subtitle = "Comparison across Control, LC with COPD, and LC without COPD",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5)
    )
  )

# Save main panel
ggsave(file.path(output_dir, "alpha_diversity_main.pdf"), 
       main_panel, width = 10, height = 10, dpi = 300)
ggsave(file.path(output_dir, "alpha_diversity_main.png"), 
       main_panel, width = 10, height = 10, dpi = 300)

# Full supplementary figure (all metrics)
if (length(alpha_plots) >= 8) {
  supp_panel <- wrap_plots(alpha_plots[1:8], ncol = 4) +
    plot_annotation(
      title = "Alpha Diversity - All Metrics",
      theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
    )
  
  ggsave(file.path(output_dir, "alpha_diversity_supplementary.pdf"), 
         supp_panel, width = 16, height = 8, dpi = 300)
  ggsave(file.path(output_dir, "alpha_diversity_supplementary.png"), 
         supp_panel, width = 16, height = 8, dpi = 300)
}

#-------------------------------------------------------------------------------
# VIOLIN PLOTS (ALTERNATIVE VISUALIZATION)
#-------------------------------------------------------------------------------
cat("Creating violin plots...\n")

plot_alpha_violin <- function(data, metric, ylab = NULL) {
  if (is.null(ylab)) ylab <- metric
  
  ggplot(data, aes(x = group, y = .data[[metric]], fill = group)) +
    geom_violin(alpha = 0.7, scale = "width") +
    geom_boxplot(width = 0.15, alpha = 0.8, outlier.shape = NA) +
    geom_jitter(width = 0.1, size = 1.5, alpha = 0.5) +
    scale_fill_manual(values = group_colors) +
    labs(x = "", y = ylab) +
    theme_publication() +
    theme(legend.position = "none")
}

violin_shannon <- plot_alpha_violin(alpha_div, "Shannon", "Shannon Diversity Index")
violin_chao1 <- plot_alpha_violin(alpha_div, "Chao1", "Chao1 Richness Estimate")

violin_panel <- violin_shannon + violin_chao1 +
  plot_annotation(title = "Alpha Diversity Distribution")

ggsave(file.path(output_dir, "alpha_diversity_violin.pdf"), 
       violin_panel, width = 10, height = 5, dpi = 300)

#-------------------------------------------------------------------------------
# RAREFACTION CURVE
#-------------------------------------------------------------------------------
cat("Creating rarefaction curves...\n")

# Prepare data for rarefaction
otu_mat <- as(otu_table(ps), "matrix")
if (taxa_are_rows(ps)) otu_mat <- t(otu_mat)

# Calculate rarefaction
raremax <- min(rowSums(otu_mat))
rarecurve_data <- rarecurve(otu_mat, step = 100, sample = raremax, 
                            tidy = TRUE)

# Add group information
rarecurve_data <- rarecurve_data %>%
  left_join(alpha_div %>% select(sample_id, group), 
            by = c("Site" = "sample_id"))

# Plot
rarefaction_plot <- ggplot(rarecurve_data, aes(x = Sample, y = Species, 
                                                 color = group, group = Site)) +
  geom_line(alpha = 0.6, linewidth = 0.5) +
  geom_vline(xintercept = raremax, linetype = "dashed", color = "gray40") +
  scale_color_manual(values = group_colors) +
  labs(x = "Sequencing Depth", y = "Observed ASVs", 
       title = "Rarefaction Curves",
       color = "Group") +
  theme_publication() +
  theme(legend.position = "right")

ggsave(file.path(output_dir, "rarefaction_curves.pdf"), 
       rarefaction_plot, width = 8, height = 6, dpi = 300)
ggsave(file.path(output_dir, "rarefaction_curves.png"), 
       rarefaction_plot, width = 8, height = 6, dpi = 300)

#-------------------------------------------------------------------------------
# SAVE RESULTS
#-------------------------------------------------------------------------------
cat("\nSaving results...\n")

# Save alpha diversity table
write.csv(alpha_div, file.path(tables_dir, "alpha_diversity_values.csv"), 
          row.names = FALSE)

# Save statistical results
write.csv(kw_results, file.path(tables_dir, "alpha_diversity_kruskal_wallis.csv"), 
          row.names = FALSE)

if (exists("pairwise_df")) {
  write.csv(pairwise_df, file.path(tables_dir, "alpha_diversity_pairwise.csv"), 
            row.names = FALSE)
}

# Save for downstream analysis
saveRDS(alpha_div, file.path(data_dir, "alpha_diversity.rds"))

cat("\n============================================\n")
cat("Alpha Diversity Analysis Complete!\n")
cat("============================================\n")
cat("\nOutput files:\n")
cat("  Figures:", output_dir, "\n")
cat("  Tables:", tables_dir, "\n")
