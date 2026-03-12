#===============================================================================
# 03_BETA_DIVERSITY.R - Beta Diversity Analysis
#===============================================================================
# Ordination, distance-based statistics, and visualization
#===============================================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(vegan)
  library(ape)
  library(ggplot2)
  library(ggpubr)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})

source("D:/LC_COPD_microbiome/R_analysis/functions/plotting_theme.R")

# Configuration
project_dir <- "D:/LC_COPD_microbiome"
data_dir <- file.path(project_dir, "R_analysis/data")
output_dir <- file.path(project_dir, "figures/beta_diversity")
tables_dir <- file.path(project_dir, "results/tables")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("============================================\n")
cat("Beta Diversity Analysis\n")
cat("============================================\n\n")

#-------------------------------------------------------------------------------
# LOAD DATA
#-------------------------------------------------------------------------------
cat("Loading data...\n")
ps <- readRDS(file.path(data_dir, "phyloseq_silva_filtered.rds"))

# Resolve polytomies: force the tree to be strictly dichotomous (binary)
# so that UniFrac edge-matrix operations do not produce length-mismatch warnings
if (!is.null(phy_tree(ps, errorIfNULL = FALSE))) {
  cat("  Resolving polytomies in the phylogenetic tree with multi2di()...\n")
  phy_tree(ps) <- multi2di(phy_tree(ps))
}

# Transform to relative abundance for some analyses
ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))

# Strip the S4 phyloseq wrapper so vegan receives a plain base data.frame
samp_data <- data.frame(sample_data(ps), stringsAsFactors = FALSE)

#-------------------------------------------------------------------------------
# CALCULATE DISTANCE MATRICES
#-------------------------------------------------------------------------------
cat("\nCalculating distance matrices...\n")

distances <- list(
  bray_curtis = phyloseq::distance(ps_rel, method = "bray"),
  jaccard = phyloseq::distance(ps_rel, method = "jaccard", binary = TRUE),
  weighted_unifrac = phyloseq::distance(ps, method = "wunifrac"),
  unweighted_unifrac = phyloseq::distance(ps, method = "unifrac")
)

cat("  Calculated: ", paste(names(distances), collapse = ", "), "\n")

#-------------------------------------------------------------------------------
# STATISTICAL TESTS
#-------------------------------------------------------------------------------
cat("\nRunning statistical tests...\n")

# PERMANOVA
permanova_results <- lapply(names(distances), function(d) {
  set.seed(42)
  test <- adonis2(distances[[d]] ~ group, data = samp_data, 
                  permutations = 999, method = d)
  data.frame(
    Distance = d,
    Df = test$Df[1],
    SumOfSqs = test$SumOfSqs[1],
    R2 = test$R2[1],
    F_value = test$F[1],
    p_value = test$`Pr(>F)`[1],
    stringsAsFactors = FALSE
  )
}) %>% bind_rows()

permanova_results$p_adjusted <- p.adjust(permanova_results$p_value, method = "BH")

cat("\nPERMANOVA Results:\n")
print(permanova_results)

# Pairwise PERMANOVA
pairwise_permanova <- function(dist_matrix, metadata, group_col) {
  groups <- unique(metadata[[group_col]])
  pairs <- combn(groups, 2, simplify = FALSE)
  
  # Build the formula as a string to avoid NSE issues with get() inside adonis2
  frm <- as.formula(paste("sub_dist ~", group_col))
  
  results <- lapply(pairs, function(pair) {
    idx <- metadata[[group_col]] %in% pair
    sub_dist <- as.dist(as.matrix(dist_matrix)[idx, idx])
    # Ensure sub_meta is a plain data.frame (strip any residual S4 class)
    sub_meta <- data.frame(metadata[idx, , drop = FALSE], stringsAsFactors = FALSE)
    
    set.seed(42)
    test <- adonis2(frm, data = sub_meta, permutations = 999)
    
    data.frame(
      Group1 = pair[1],
      Group2 = pair[2],
      R2 = test$R2[1],
      F_value = test$F[1],
      p_value = test$`Pr(>F)`[1]
    )
  }) %>% bind_rows()
  
  results$p_adjusted <- p.adjust(results$p_value, method = "BH")
  return(results)
}

pairwise_results <- lapply(names(distances), function(d) {
  res <- pairwise_permanova(distances[[d]], samp_data, "group")
  res$Distance <- d
  return(res)
}) %>% bind_rows()

cat("\nPairwise PERMANOVA:\n")
print(pairwise_results)

# Beta dispersion test
betadisper_results <- lapply(names(distances), function(d) {
  bd <- betadisper(distances[[d]], samp_data$group)
  test <- permutest(bd, permutations = 999)
  
  data.frame(
    Distance = d,
    F_value = test$statistic,
    p_value = test$tab$`Pr(>F)`[1]
  )
}) %>% bind_rows()

cat("\nBeta Dispersion Tests:\n")
print(betadisper_results)

#-------------------------------------------------------------------------------
# PCoA ORDINATIONS
#-------------------------------------------------------------------------------
cat("\nPerforming PCoA ordinations...\n")

pcoa_plots <- list()

for (d in names(distances)) {
  # Calculate PCoA
  pcoa <- ordinate(ps, method = "PCoA", distance = distances[[d]])
  
  # Get variance explained
  var_explained <- round(pcoa$values$Relative_eig[1:2] * 100, 1)
  
  # Create plot
  p <- plot_ordination(ps, pcoa, color = "group") +
    geom_point(size = 3.5, alpha = 0.85) +
    stat_ellipse(aes(color = group), type = "norm", level = 0.95, 
                 linewidth = 0.8, linetype = "dashed") +
    scale_color_manual(values = group_colors, name = "Group") +
    labs(
      x = sprintf("PCoA1 (%.1f%%)", var_explained[1]),
      y = sprintf("PCoA2 (%.1f%%)", var_explained[2]),
      title = gsub("_", "-", tools::toTitleCase(d))
    ) +
    theme_publication() +
    theme(legend.position = "right")
  
  # Add PERMANOVA annotation
  perm_p <- permanova_results$p_value[permanova_results$Distance == d]
  perm_r2 <- permanova_results$R2[permanova_results$Distance == d]
  p <- p + annotate("text", x = Inf, y = Inf, 
                    label = sprintf("PERMANOVA\nR² = %.3f\np = %.3f", perm_r2, perm_p),
                    hjust = 1.1, vjust = 1.5, size = 3, fontface = "italic")
  
  pcoa_plots[[d]] <- p
  
  # Save individual plot
  ggsave(file.path(output_dir, paste0("pcoa_", d, ".pdf")), 
         p, width = 8, height = 6, dpi = 300)
}

#-------------------------------------------------------------------------------
# COMBINED PCoA FIGURE
#-------------------------------------------------------------------------------
cat("Creating combined PCoA figure...\n")

combined_pcoa <- (pcoa_plots$bray_curtis + pcoa_plots$jaccard) /
                 (pcoa_plots$weighted_unifrac + pcoa_plots$unweighted_unifrac) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Beta Diversity Ordinations",
    theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
  ) &
  theme(legend.position = "bottom")

ggsave(file.path(output_dir, "beta_diversity_pcoa_combined.pdf"), 
       combined_pcoa, width = 12, height = 10, dpi = 300)
ggsave(file.path(output_dir, "beta_diversity_pcoa_combined.png"), 
       combined_pcoa, width = 12, height = 10, dpi = 300)

#-------------------------------------------------------------------------------
# NMDS ORDINATION
#-------------------------------------------------------------------------------
cat("Performing NMDS ordination...\n")

set.seed(42)
nmds <- ordinate(ps_rel, method = "NMDS", distance = "bray", trymax = 100)

nmds_plot <- plot_ordination(ps_rel, nmds, color = "group") +
  geom_point(size = 3.5, alpha = 0.85) +
  stat_ellipse(aes(color = group), type = "norm", level = 0.95,
               linewidth = 0.8, linetype = "dashed") +
  scale_color_manual(values = group_colors) +
  labs(title = "NMDS (Bray-Curtis)",
       subtitle = sprintf("Stress = %.3f", nmds$stress)) +
  theme_publication() +
  theme(legend.position = "right")

ggsave(file.path(output_dir, "nmds_bray_curtis.pdf"), 
       nmds_plot, width = 8, height = 6, dpi = 300)

#-------------------------------------------------------------------------------
# DISTANCE BOXPLOTS (WITHIN VS BETWEEN GROUPS)
#-------------------------------------------------------------------------------
cat("Creating distance comparison boxplots...\n")

# Function to categorize distances
categorize_distances <- function(dist_matrix, metadata, group_col) {
  dist_mat <- as.matrix(dist_matrix)
  groups <- metadata[[group_col]]
  n <- nrow(dist_mat)
  
  results <- data.frame()
  
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      g1 <- groups[i]
      g2 <- groups[j]
      
      if (g1 == g2) {
        category <- paste0("Within ", g1)
      } else {
        category <- paste(sort(c(as.character(g1), as.character(g2))), collapse = " vs ")
      }
      
      results <- rbind(results, data.frame(
        Distance = dist_mat[i, j],
        Category = category,
        Type = ifelse(g1 == g2, "Within", "Between")
      ))
    }
  }
  
  return(results)
}

# Create boxplot for Bray-Curtis
dist_categories <- categorize_distances(distances$bray_curtis, samp_data, "group")

dist_boxplot <- ggplot(dist_categories, aes(x = Category, y = Distance, fill = Type)) +
  geom_boxplot(alpha = 0.7, outlier.size = 1) +
  scale_fill_manual(values = c("Within" = "#4477AA", "Between" = "#EE6677")) +
  labs(x = "", y = "Bray-Curtis Distance",
       title = "Distance Comparison: Within vs Between Groups") +
  theme_publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(output_dir, "distance_comparison_boxplot.pdf"), 
       dist_boxplot, width = 10, height = 6, dpi = 300)

#-------------------------------------------------------------------------------
# SAVE RESULTS
#-------------------------------------------------------------------------------
cat("\nSaving results...\n")

# Save statistical results
write.csv(permanova_results, 
          file.path(tables_dir, "beta_diversity_permanova.csv"), 
          row.names = FALSE)

write.csv(pairwise_results, 
          file.path(tables_dir, "beta_diversity_pairwise_permanova.csv"), 
          row.names = FALSE)

write.csv(betadisper_results, 
          file.path(tables_dir, "beta_diversity_dispersion.csv"), 
          row.names = FALSE)

# Save distance matrices
for (d in names(distances)) {
  write.csv(as.matrix(distances[[d]]), 
            file.path(tables_dir, paste0("distance_matrix_", d, ".csv")))
}

cat("\n============================================\n")
cat("Beta Diversity Analysis Complete!\n")
cat("============================================\n")
