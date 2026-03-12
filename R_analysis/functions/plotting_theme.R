#===============================================================================
# PLOTTING_THEME.R - Custom Plotting Functions and Theme
#===============================================================================
# Consistent styling for publication-quality microbiome plots
#===============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
  library(viridis)
  library(RColorBrewer)
})

#-------------------------------------------------------------------------------
# COLOR PALETTES
#-------------------------------------------------------------------------------

# Main group colors (colorblind-friendly)
group_colors <- c(
  "Control" = "#4477AA",      # Blue
  "LC_COPD" = "#EE6677",      # Red
  "LC_woCOPD" = "#228833"     # Green
)

# Alternative group colors
group_colors_alt <- c(
  "Control" = "#0077BB",
  "LC_COPD" = "#CC3311", 
  "LC_woCOPD" = "#009988"
)

# Extended palette for taxa (12 colors)
taxa_colors <- c(
  "#E64B35", "#4DBBD5", "#00A087", "#3C5488",
  "#F39B7F", "#8491B4", "#91D1C2", "#DC0000",
  "#7E6148", "#B09C85", "#00468B", "#42B540"
)

# Large palette for many taxa (20+ colors)
taxa_colors_extended <- c(
  colorRampPalette(brewer.pal(12, "Paired"))(20)
)

# Sequential palette for heatmaps
heatmap_colors <- viridis::viridis(100)
heatmap_colors_alt <- colorRampPalette(c("white", "#3182BD", "#08519C"))(100)

# Diverging palette
diverging_colors <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100)

#-------------------------------------------------------------------------------
# GGPLOT2 PUBLICATION THEME
#-------------------------------------------------------------------------------

theme_publication <- function(base_size = 12, base_family = "") {
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      # Panel
      panel.background = element_rect(fill = "white", color = NA),
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
      panel.grid.major = element_line(color = "gray90", linewidth = 0.25),
      panel.grid.minor = element_blank(),
      
      # Axes
      axis.title = element_text(size = rel(1.0), face = "bold"),
      axis.text = element_text(size = rel(0.9), color = "black"),
      axis.ticks = element_line(color = "black", linewidth = 0.25),
      axis.line = element_blank(),
      
      # Legend
      legend.background = element_rect(fill = "white", color = NA),
      legend.key = element_rect(fill = "white", color = NA),
      legend.key.size = unit(0.8, "lines"),
      legend.title = element_text(size = rel(0.9), face = "bold"),
      legend.text = element_text(size = rel(0.85)),
      legend.position = "right",
      
      # Strips (for facets)
      strip.background = element_rect(fill = "gray95", color = "black", linewidth = 0.5),
      strip.text = element_text(size = rel(0.9), face = "bold"),
      
      # Plot title
      plot.title = element_text(size = rel(1.2), face = "bold", hjust = 0.5,
                                margin = margin(b = 10)),
      plot.subtitle = element_text(size = rel(1.0), hjust = 0.5,
                                   margin = margin(b = 10)),
      plot.caption = element_text(size = rel(0.8), hjust = 1),
      
      # Margins
      plot.margin = margin(10, 10, 10, 10)
    )
}

# Minimal theme for clean plots
theme_minimal_pub <- function(base_size = 12) {
  theme_minimal(base_size = base_size) %+replace%
    theme(
      panel.grid.major = element_line(color = "gray90", linewidth = 0.25),
      panel.grid.minor = element_blank(),
      axis.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
}

#-------------------------------------------------------------------------------
# HELPER FUNCTIONS
#-------------------------------------------------------------------------------

# Save publication-quality figure
save_publication_plot <- function(plot, filename, width = 8, height = 6, 
                                   dpi = 300, formats = c("pdf", "png")) {
  for (fmt in formats) {
    outfile <- paste0(tools::file_path_sans_ext(filename), ".", fmt)
    ggsave(outfile, plot, width = width, height = height, dpi = dpi)
    cat(sprintf("Saved: %s\n", outfile))
  }
}

# Format p-values for display
format_pvalue <- function(p) {
  case_when(
    p < 0.001 ~ "< 0.001",
    p < 0.01 ~ sprintf("%.3f", p),
    p < 0.05 ~ sprintf("%.3f", p),
    TRUE ~ sprintf("%.2f", p)
  )
}

# Add significance stars
add_significance <- function(p) {
  case_when(
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ "ns"
  )
}

# Get top N taxa names
get_top_taxa <- function(ps, n = 10, rank = "Genus") {
  tax <- tax_table(ps)
  otu <- otu_table(ps)
  
  if (taxa_are_rows(ps)) {
    taxa_sums <- rowSums(otu)
  } else {
    taxa_sums <- colSums(otu)
  }
  
  top_idx <- order(taxa_sums, decreasing = TRUE)[1:n]
  top_names <- tax[top_idx, rank]
  
  return(top_names)
}

# Aggregate taxa by rank
aggregate_taxa <- function(ps, rank = "Genus") {
  ps_glom <- tax_glom(ps, taxrank = rank)
  return(ps_glom)
}

# Calculate relative abundance
to_relative_abundance <- function(ps) {
  transform_sample_counts(ps, function(x) x / sum(x) * 100)
}

#-------------------------------------------------------------------------------
# ORDINATION PLOT FUNCTION
#-------------------------------------------------------------------------------

plot_ordination_custom <- function(ps, method = "PCoA", distance = "bray",
                                    color_by = "group", shape_by = NULL,
                                    ellipse = TRUE, title = NULL) {
  
  # Calculate ordination
  ord <- ordinate(ps, method = method, distance = distance)
  
  # Create plot
  p <- plot_ordination(ps, ord, color = color_by, shape = shape_by) +
    geom_point(size = 3, alpha = 0.8) +
    scale_color_manual(values = group_colors) +
    theme_publication()
  
  # Add ellipses
  if (ellipse) {
    p <- p + stat_ellipse(aes(color = .data[[color_by]]), 
                          type = "norm", level = 0.95, 
                          linewidth = 0.8, linetype = "dashed")
  }
  
  # Add title
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  
  return(p)
}

#-------------------------------------------------------------------------------
# BARPLOT FUNCTION
#-------------------------------------------------------------------------------

plot_taxa_barplot <- function(ps, rank = "Phylum", n_taxa = 10, 
                               group_by = "group", title = NULL) {
  
  # Aggregate to rank
  ps_rank <- tax_glom(ps, taxrank = rank)
  
  # Transform to relative abundance
  ps_rel <- transform_sample_counts(ps_rank, function(x) x / sum(x) * 100)
  
  # Melt to long format
  df <- psmelt(ps_rel)
  
  # Get top taxa
  top_taxa <- df %>%
    group_by(across(all_of(rank))) %>%
    summarise(mean_abund = mean(Abundance), .groups = "drop") %>%
    arrange(desc(mean_abund)) %>%
    slice_head(n = n_taxa) %>%
    pull(all_of(rank))
  
  # Label others
  df <- df %>%
    mutate(Taxon = ifelse(get(rank) %in% top_taxa, get(rank), "Other")) %>%
    mutate(Taxon = factor(Taxon, levels = c(top_taxa, "Other")))
  
  # Aggregate by taxon
  df_agg <- df %>%
    group_by(Sample, Taxon, across(all_of(group_by))) %>%
    summarise(Abundance = sum(Abundance), .groups = "drop")
  
  # Plot
  p <- ggplot(df_agg, aes(x = Sample, y = Abundance, fill = Taxon)) +
    geom_bar(stat = "identity", width = 0.9) +
    facet_grid(~ get(group_by), scales = "free_x", space = "free_x") +
    scale_fill_manual(values = c(taxa_colors[1:n_taxa], "gray70")) +
    labs(x = "", y = "Relative Abundance (%)", fill = rank) +
    theme_publication() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      legend.position = "right"
    )
  
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  
  return(p)
}

cat("Plotting functions loaded successfully.\n")
cat("Available color palettes: group_colors, taxa_colors, heatmap_colors\n")
cat("Available themes: theme_publication(), theme_minimal_pub()\n")
