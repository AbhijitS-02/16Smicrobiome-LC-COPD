# ===============================================================================
# 06_PICRUST_VISUALIZATION_V2.R - PICRUSt2 Functional Analysis Visualization
# ===============================================================================
# Optimized version of 06_picrust_visualization.R
# Changes from v1:
#   1. Vectorized Kruskal-Wallis testing (no rbind loop — O(1) memory growth)
#   2. Removed orphaned KO data loading (ko_data was never used)
#   3. Fixed barplot sorting (by q-value instead of mean abundance across groups)
#   4. Safe heatmap rendering via tryCatch (prevents orphan PDF devices)
# ===============================================================================

suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(patchwork)
    library(pheatmap)
    library(viridis)
    library(vegan)
})

source("D:/LC_COPD_microbiome/R_analysis/functions/plotting_theme.R")

# Configuration
project_dir <- "D:/LC_COPD_microbiome"
data_dir <- file.path(project_dir, "R_analysis/data")
picrust_dir <- file.path(project_dir, "results/picrust2")
output_dir <- file.path(project_dir, "figures/functional")
tables_dir <- file.path(project_dir, "results/tables")
metadata_file <- file.path(project_dir, "metadata/sample_metadata.tsv")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("============================================\n")
cat("PICRUSt2 Functional Analysis (v2)\n")
cat("============================================\n\n")

#-------------------------------------------------------------------------------
# LOAD DATA
#-------------------------------------------------------------------------------
cat("Loading PICRUSt2 data...\n")

# Load metadata
metadata <- read.table(metadata_file,
    header = TRUE, sep = "\t",
    row.names = 1, stringsAsFactors = FALSE
)
metadata$group <- factor(metadata$group, levels = c("Control", "LC_COPD", "LC_woCOPD"))

# Load pathway predictions
pathway_file <- file.path(tables_dir, "pathways_with_descriptions.tsv")

if (!file.exists(pathway_file)) {
    # Try alternative locations
    pathway_file <- file.path(picrust_dir, "pathways_with_descriptions.tsv")
}

if (!file.exists(pathway_file)) {
    pathway_file <- file.path(picrust_dir, "picrust2_out/pathways_out/path_abun_unstrat.tsv")
    if (file.exists(paste0(pathway_file, ".gz"))) {
        system(paste("gunzip -k", paste0(pathway_file, ".gz")))
    }
}

if (file.exists(pathway_file)) {
    pathways <- read.table(pathway_file,
        header = TRUE, sep = "\t",
        row.names = 1, check.names = FALSE,
        quote = "", comment.char = ""
    )
    cat(sprintf(
        "  Loaded %d pathways across %d samples\n",
        nrow(pathways), ncol(pathways)
    ))
} else {
    stop("Pathway file not found. Please run PICRUSt2 first.")
}

# NOTE: KO data loading removed — ko_data was never referenced downstream.
#       If KO-level analysis is needed in the future, replicate the pathway
#       pipeline (PCA, KW testing, heatmap) for the ko_data matrix.

#-------------------------------------------------------------------------------
# PROCESS PATHWAY DATA
#-------------------------------------------------------------------------------
cat("\nProcessing pathway data...\n")

# Separate description column if present
if ("description" %in% colnames(pathways)) {
    pathway_desc <- pathways$description
    names(pathway_desc) <- rownames(pathways)
    pathways <- pathways[, !colnames(pathways) %in% "description"]
}

# Ensure samples match metadata
common_samples <- intersect(colnames(pathways), rownames(metadata))

cat(sprintf(
    "  Pathway columns (first 5): %s\n",
    paste(head(colnames(pathways), 5), collapse = ", ")
))
cat(sprintf(
    "  Metadata rows (first 5):   %s\n",
    paste(head(rownames(metadata), 5), collapse = ", ")
))
cat(sprintf("  Matching samples: %d\n", length(common_samples)))

if (length(common_samples) == 0) {
    stop(
        "No matching sample IDs between pathway data columns and metadata rows!\n",
        "  Check that sample IDs are consistent between files."
    )
}

pathways <- pathways[, common_samples]
metadata <- metadata[common_samples, , drop = FALSE]

# Convert to relative abundance
pathway_rel <- apply(pathways, 2, function(x) x / sum(x) * 100)

cat(sprintf("  Final: %d pathways, %d samples\n", nrow(pathway_rel), ncol(pathway_rel)))

#-------------------------------------------------------------------------------
# PATHWAY PCA
#-------------------------------------------------------------------------------
cat("\nPerforming PCA on pathway abundances...\n")

# Transpose for PCA
pathway_t <- t(pathway_rel)

# PCA
pca_result <- prcomp(pathway_t, scale. = TRUE, center = TRUE)

# Variance explained
var_explained <- round(pca_result$sdev^2 / sum(pca_result$sdev^2) * 100, 1)

# Create plot data
pca_df <- data.frame(
    PC1 = pca_result$x[, 1],
    PC2 = pca_result$x[, 2],
    group = metadata$group
)

# PCA plot
pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = group)) +
    geom_point(size = 3.5, alpha = 0.85) +
    stat_ellipse(type = "norm", level = 0.95, linewidth = 0.8, linetype = "dashed") +
    scale_color_manual(values = group_colors) +
    labs(
        x = sprintf("PC1 (%.1f%%)", var_explained[1]),
        y = sprintf("PC2 (%.1f%%)", var_explained[2]),
        title = "PCA of Predicted Pathway Abundances",
        color = "Group"
    ) +
    theme_publication()

ggsave(file.path(output_dir, "pathway_pca.pdf"),
    pca_plot,
    width = 8, height = 6, dpi = 300
)
ggsave(file.path(output_dir, "pathway_pca.png"),
    pca_plot,
    width = 8, height = 6, dpi = 300
)

#-------------------------------------------------------------------------------
# PERMANOVA ON PATHWAYS
#-------------------------------------------------------------------------------
cat("Running PERMANOVA on pathway abundances...\n")

# Calculate Bray-Curtis distance
pathway_dist <- vegdist(pathway_t, method = "bray")

# PERMANOVA
set.seed(42)
permanova <- adonis2(pathway_dist ~ group, data = metadata, permutations = 200)

cat("  PERMANOVA Results:\n")
print(permanova)

# Add to PCA plot
pca_plot_annotated <- pca_plot +
    annotate("text",
        x = Inf, y = Inf,
        label = sprintf(
            "PERMANOVA\nR² = %.3f\np = %.3f",
            permanova$R2[1], permanova$`Pr(>F)`[1]
        ),
        hjust = 1.1, vjust = 1.1, size = 3.5, fontface = "italic"
    )

ggsave(file.path(output_dir, "pathway_pca_annotated.pdf"),
    pca_plot_annotated,
    width = 8, height = 6, dpi = 300
)

#-------------------------------------------------------------------------------
# TOP DIFFERENTIAL PATHWAYS  (vectorized — no rbind loop)
#-------------------------------------------------------------------------------
cat("\nIdentifying differentially abundant pathways...\n")

# --- FIX 1: Vectorized Kruskal-Wallis testing ---
# Using apply() over rows avoids the O(n²) cost of rbind-in-a-loop.
# Each iteration returns a single p-value; the data.frame is built once.

p_values <- apply(pathway_rel, 1, function(x) {
    kruskal.test(as.numeric(x) ~ metadata$group)$p.value
})

# Calculate means per group efficiently using tapply on column indices
group_idx <- split(seq_along(metadata$group), metadata$group)
means_control <- rowMeans(pathway_rel[, group_idx[["Control"]], drop = FALSE])
means_lc_copd <- rowMeans(pathway_rel[, group_idx[["LC_COPD"]], drop = FALSE])
means_lc_wocopd <- rowMeans(pathway_rel[, group_idx[["LC_woCOPD"]], drop = FALSE])

# Build results data.frame in one shot
kw_results <- data.frame(
    Pathway = rownames(pathway_rel),
    Description = if (exists("pathway_desc")) {
        pathway_desc[rownames(pathway_rel)]
    } else {
        rownames(pathway_rel)
    },
    Mean_Control = means_control,
    Mean_LC_COPD = means_lc_copd,
    Mean_LC_woCOPD = means_lc_wocopd,
    p_value = p_values,
    stringsAsFactors = FALSE
)

# Adjust p-values
kw_results$q_value <- p.adjust(kw_results$p_value, method = "BH")
kw_results <- kw_results %>% arrange(p_value)

# Significant pathways (using raw p-value)
sig_pathways <- kw_results %>% filter(p_value < 0.05)
cat(sprintf("  Significant pathways (p < 0.05): %d\n", nrow(sig_pathways)))

#-------------------------------------------------------------------------------
# HEATMAP OF TOP 30 SIGNIFICANT PATHWAYS  (p < 0.05)
#-------------------------------------------------------------------------------
cat("Creating pathway heatmap...\n")

# Use top 30 significant pathways (already sorted by p-value in sig_pathways)
if (nrow(sig_pathways) > 0) {
    top30_sig <- sig_pathways %>% slice_head(n = 30)
    top30_pathways <- pathway_rel[top30_sig$Pathway, ]

    # Shorten pathway names if too long
    row_labels <- rownames(top30_pathways)
    if (exists("pathway_desc")) {
        row_labels <- pathway_desc[rownames(top30_pathways)]
    }
    row_labels <- substr(row_labels, 1, 50) # Truncate long names

    # Annotation
    annotation_col <- data.frame(Group = metadata$group, row.names = rownames(metadata))

    # Log transform for visualization
    heatmap_mat <- log10(top30_pathways + 0.001)
    rownames(heatmap_mat) <- row_labels

    # Safe heatmap rendering via tryCatch
    tryCatch(
        {
            pdf(file.path(output_dir, "pathway_heatmap_top30.pdf"), width = 14, height = 10)
            pheatmap(
                heatmap_mat,
                annotation_col = annotation_col,
                annotation_colors = list(Group = group_colors),
                color = viridis(100),
                cluster_rows = TRUE,
                cluster_cols = TRUE,
                fontsize_row = 7,
                fontsize_col = 8,
                main = sprintf(
                    "Top %d Significant MetaCyc Pathways (p < 0.05)",
                    nrow(top30_sig)
                )
            )
            dev.off()
            cat(sprintf("  Heatmap saved (%d pathways).\n", nrow(top30_sig)))
        },
        error = function(e) {
            cat(sprintf("  WARNING: Heatmap generation failed (%s). Skipping.\n", e$message))
            try(dev.off(), silent = TRUE)
        }
    )
} else {
    cat("  No significant pathways for heatmap. Skipping.\n")
}

#-------------------------------------------------------------------------------
# BARPLOT OF SIGNIFICANT PATHWAYS  (with pairwise significance annotations)
#-------------------------------------------------------------------------------
if (nrow(sig_pathways) > 0) {
    cat("Creating significant pathway barplot...\n")

    # Take top 15 significant pathways (sorted by p-value)
    top_sig <- sig_pathways %>% slice_head(n = 15)

    # Pivot to long format for plotting
    plot_data <- top_sig %>%
        mutate(Pathway_Short = substr(Pathway, 1, 40)) %>%
        mutate(Pathway_Short = factor(Pathway_Short, levels = rev(Pathway_Short))) %>%
        pivot_longer(
            cols = starts_with("Mean_"),
            names_to = "Group", values_to = "Abundance"
        ) %>%
        mutate(
            Group = gsub("Mean_", "", Group),
            Group = factor(Group, levels = c("Control", "LC_COPD", "LC_woCOPD"))
        )

    # --- Compute pairwise Wilcoxon tests for significance brackets ---
    comparisons <- list(
        c("Control", "LC_COPD"),
        c("Control", "LC_woCOPD"),
        c("LC_COPD", "LC_woCOPD")
    )
    # Symbols for each comparison
    comp_symbols <- c("*", "^", "#")

    # For each pathway, run pairwise Wilcoxon on the raw relative abundances
    annot_list <- list()
    for (i in seq_len(nrow(top_sig))) {
        pw <- top_sig$Pathway[i]
        pw_short <- substr(pw, 1, 40)
        pw_vals <- as.numeric(pathway_rel[pw, ])
        max_abund <- max(
            top_sig$Mean_Control[i], top_sig$Mean_LC_COPD[i],
            top_sig$Mean_LC_woCOPD[i]
        )

        for (j in seq_along(comparisons)) {
            g1 <- comparisons[[j]][1]
            g2 <- comparisons[[j]][2]
            vals1 <- pw_vals[metadata$group == g1]
            vals2 <- pw_vals[metadata$group == g2]

            p_val <- tryCatch(
                wilcox.test(vals1, vals2)$p.value,
                error = function(e) NA
            )

            # Determine symbol
            if (is.na(p_val)) {
                sym <- "ns"
            } else if (p_val < 0.001) {
                sym <- paste0(comp_symbols[j], comp_symbols[j], comp_symbols[j])
            } else if (p_val < 0.01) {
                sym <- paste0(comp_symbols[j], comp_symbols[j])
            } else if (p_val < 0.05) {
                sym <- comp_symbols[j]
            } else {
                sym <- "ns"
            }

            annot_list[[length(annot_list) + 1]] <- data.frame(
                Pathway_Short = pw_short,
                Comparison = paste0(g1, " vs ", g2),
                comp_idx = j,
                symbol = sym,
                p_value = p_val,
                max_abund = max_abund,
                stringsAsFactors = FALSE
            )
        }
    }
    annot_df <- do.call(rbind, annot_list)
    annot_df$Pathway_Short <- factor(annot_df$Pathway_Short,
        levels = levels(plot_data$Pathway_Short)
    )

    # Position symbols using a GLOBAL offset so low-abundance pathways
    # don't have their symbols piled on top of each other.
    global_max <- max(annot_df$max_abund, na.rm = TRUE)
    step_size <- global_max * 0.08 # fixed gap between comparison columns

    annot_df <- annot_df %>%
        mutate(
            # Place each symbol at: pathway's max bar + fixed offset per comparison
            y_pos = max_abund + step_size * comp_idx
        )

    # Build the base barplot
    pathway_bar <- ggplot(plot_data, aes(
        x = Pathway_Short,
        y = Abundance, fill = Group
    )) +
        geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
        coord_flip() +
        scale_fill_manual(values = group_colors) +
        labs(
            x = "", y = "Mean Relative Abundance (%)",
            title = "Top Differentially Abundant Pathways (p < 0.05)"
        ) +
        theme_publication() +
        theme(axis.text.y = element_text(size = 8))

    # Add significance annotations as text labels
    # For a horizontal barplot (coord_flip), x-axis is the pathway, y-axis is abundance
    # We place the symbols to the right of the tallest bar, staggered per comparison
    pathway_bar <- pathway_bar +
        geom_text(
            data = annot_df,
            aes(
                x = Pathway_Short,
                y = y_pos,
                label = symbol
            ),
            inherit.aes = FALSE,
            size = 3.5, fontface = "bold", hjust = 0.5, vjust = 0.5
        )

    # Expand x-axis to make room for annotations
    pathway_bar <- pathway_bar +
        scale_y_continuous(expand = expansion(mult = c(0.02, 0.35)))

    # Add a legend for the significance symbols at the bottom
    sig_legend <- paste0(
        "Pairwise Wilcoxon: * = Control vs LC_COPD, ",
        "^ = Control vs LC_woCOPD, ",
        "# = LC_COPD vs LC_woCOPD\n",
        "Significance: one symbol (p < 0.05), ",
        "two symbols (p < 0.01), ",
        "three symbols (p < 0.001), ",
        "ns = not significant (p >= 0.05)"
    )
    pathway_bar <- pathway_bar +
        labs(caption = sig_legend) +
        theme(plot.caption = element_text(hjust = 0, size = 7, face = "italic"))

    ggsave(file.path(output_dir, "significant_pathways_barplot.pdf"),
        pathway_bar,
        width = 14, height = 8, dpi = 300
    )
    ggsave(file.path(output_dir, "significant_pathways_barplot.png"),
        pathway_bar,
        width = 14, height = 8, dpi = 300
    )
    cat("  Barplot with pairwise annotations saved.\n")
}

#-------------------------------------------------------------------------------
# SAVE RESULTS
#-------------------------------------------------------------------------------
cat("\nSaving results...\n")

write.csv(kw_results, file.path(tables_dir, "pathway_differential_analysis.csv"),
    row.names = FALSE
)

if (nrow(sig_pathways) > 0) {
    write.csv(sig_pathways, file.path(tables_dir, "pathway_significant.csv"),
        row.names = FALSE
    )
}

# Save PERMANOVA result
permanova_df <- data.frame(
    Source = rownames(permanova),
    permanova
)
write.csv(permanova_df, file.path(tables_dir, "pathway_permanova.csv"),
    row.names = FALSE
)

cat("\n============================================\n")
cat("PICRUSt2 Visualization Complete! (v2)\n")
cat("============================================\n")
