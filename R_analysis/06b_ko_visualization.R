# ===============================================================================
# 06b_KO_VISUALIZATION.R - PICRUSt2 KO (KEGG Orthology) Analysis Visualization
# ===============================================================================
# Mirrors the pathway visualization pipeline (06_picrust_visualization.R)
# but operates on KEGG Orthology (KO) metagenome predictions.
#
# Outputs:
#   - PCA of KO abundances (with PERMANOVA annotation)
#   - Heatmap of top 30 significant KOs (p < 0.05)
#   - Barplot of top 15 significant KOs with pairwise Wilcoxon annotations
#   - CSV tables of differential analysis results
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
cat("PICRUSt2 KO (KEGG Orthology) Analysis\n")
cat("============================================\n\n")

#-------------------------------------------------------------------------------
# LOAD DATA
#-------------------------------------------------------------------------------
cat("Loading KO data...\n")

# Load metadata
metadata <- read.table(metadata_file,
    header = TRUE, sep = "\t",
    row.names = 1, stringsAsFactors = FALSE
)
metadata$group <- factor(metadata$group, levels = c("Control", "LC_COPD", "LC_woCOPD"))

# Load KO metagenome predictions
ko_file <- file.path(tables_dir, "KO_metagenome_with_descriptions.tsv")

if (!file.exists(ko_file)) {
    ko_file <- file.path(picrust_dir, "KO_metagenome_with_descriptions.tsv")
}

if (!file.exists(ko_file)) {
    ko_file <- file.path(picrust_dir, "picrust2_out/KO_metagenome_out/pred_metagenome_unstrat.tsv")
    if (file.exists(paste0(ko_file, ".gz"))) {
        system(paste("gunzip -k", paste0(ko_file, ".gz")))
    }
}

if (file.exists(ko_file)) {
    ko_data <- read.table(ko_file,
        header = TRUE, sep = "\t",
        row.names = 1, check.names = FALSE,
        quote = "", comment.char = ""
    )
    cat(sprintf(
        "  Loaded %d KOs across %d columns\n",
        nrow(ko_data), ncol(ko_data)
    ))
} else {
    stop("KO file not found. Please run PICRUSt2 first.")
}

#-------------------------------------------------------------------------------
# PROCESS KO DATA
#-------------------------------------------------------------------------------
cat("\nProcessing KO data...\n")

# Separate description column if present
if ("description" %in% colnames(ko_data)) {
    ko_desc <- ko_data$description
    names(ko_desc) <- rownames(ko_data)
    ko_data <- ko_data[, !colnames(ko_data) %in% "description"]
}

# Ensure samples match metadata
common_samples <- intersect(colnames(ko_data), rownames(metadata))

cat(sprintf(
    "  KO columns (first 5): %s\n",
    paste(head(colnames(ko_data), 5), collapse = ", ")
))
cat(sprintf(
    "  Metadata rows (first 5):   %s\n",
    paste(head(rownames(metadata), 5), collapse = ", ")
))
cat(sprintf("  Matching samples: %d\n", length(common_samples)))

if (length(common_samples) == 0) {
    stop(
        "No matching sample IDs between KO data columns and metadata rows!\n",
        "  Check that sample IDs are consistent between files."
    )
}

ko_data <- ko_data[, common_samples]
metadata <- metadata[common_samples, , drop = FALSE]

# Convert to numeric (precaution)
ko_data <- as.data.frame(lapply(ko_data, as.numeric))
rownames(ko_data) <- names(ko_desc) # restore rownames after lapply

# Convert to relative abundance
ko_rel <- apply(ko_data, 2, function(x) x / sum(x) * 100)

cat(sprintf("  Final: %d KOs, %d samples\n", nrow(ko_rel), ncol(ko_rel)))

#-------------------------------------------------------------------------------
# KO PCA
#-------------------------------------------------------------------------------
cat("\nPerforming PCA on KO abundances...\n")

# Remove zero-variance rows before PCA (KOs with 0 in all samples)
ko_var_all <- apply(ko_rel, 1, var)
ko_rel_nzv <- ko_rel[ko_var_all > 0, ]
cat(sprintf("  Removed %d zero-variance KOs for PCA\n", sum(ko_var_all == 0)))

# Transpose for PCA
ko_t <- t(ko_rel_nzv)

# PCA
pca_result <- prcomp(ko_t, scale. = TRUE, center = TRUE)

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
        title = "PCA of Predicted KO Abundances",
        color = "Group"
    ) +
    theme_publication()

ggsave(file.path(output_dir, "ko_pca.pdf"),
    pca_plot,
    width = 8, height = 6, dpi = 300
)
ggsave(file.path(output_dir, "ko_pca.png"),
    pca_plot,
    width = 8, height = 6, dpi = 300
)

#-------------------------------------------------------------------------------
# PERMANOVA ON KOs
#-------------------------------------------------------------------------------
cat("Running PERMANOVA on KO abundances...\n")

# Calculate Bray-Curtis distance
ko_dist <- vegdist(ko_t, method = "bray")

# PERMANOVA
set.seed(42)
permanova <- adonis2(ko_dist ~ group, data = metadata, permutations = 200)

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

ggsave(file.path(output_dir, "ko_pca_annotated.pdf"),
    pca_plot_annotated,
    width = 8, height = 6, dpi = 300
)

#-------------------------------------------------------------------------------
# TOP DIFFERENTIAL KOs  (vectorized)
#-------------------------------------------------------------------------------
cat("\nIdentifying differentially abundant KOs...\n")

# Vectorized Kruskal-Wallis testing
p_values <- apply(ko_rel, 1, function(x) {
    kruskal.test(as.numeric(x) ~ metadata$group)$p.value
})

# Calculate means per group
group_idx <- split(seq_along(metadata$group), metadata$group)
means_control <- rowMeans(ko_rel[, group_idx[["Control"]], drop = FALSE])
means_lc_copd <- rowMeans(ko_rel[, group_idx[["LC_COPD"]], drop = FALSE])
means_lc_wocopd <- rowMeans(ko_rel[, group_idx[["LC_woCOPD"]], drop = FALSE])

# Build results data.frame
kw_results <- data.frame(
    KO = rownames(ko_rel),
    Description = if (exists("ko_desc")) {
        ko_desc[rownames(ko_rel)]
    } else {
        rownames(ko_rel)
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

# Significant KOs (using raw p-value)
sig_kos <- kw_results %>% filter(p_value < 0.05)
cat(sprintf("  Significant KOs (p < 0.05): %d\n", nrow(sig_kos)))

#-------------------------------------------------------------------------------
# HEATMAP OF TOP 30 SIGNIFICANT KOs  (p < 0.05)
#-------------------------------------------------------------------------------
cat("Creating KO heatmap...\n")

# Use top 30 significant KOs (already sorted by p-value)
if (nrow(sig_kos) > 0) {
    top30_sig <- sig_kos %>% slice_head(n = 30)
    top30_kos <- ko_rel[top30_sig$KO, ]

    # Use description as row labels if available
    row_labels <- rownames(top30_kos)
    if (exists("ko_desc")) {
        row_labels <- ko_desc[rownames(top30_kos)]
    }
    row_labels <- substr(row_labels, 1, 60) # Truncate long names

    # Annotation
    annotation_col <- data.frame(Group = metadata$group, row.names = rownames(metadata))

    # Log transform for visualization
    heatmap_mat <- log10(top30_kos + 0.001)
    rownames(heatmap_mat) <- row_labels

    # Safe heatmap rendering via tryCatch
    tryCatch(
        {
            pdf(file.path(output_dir, "ko_heatmap_top30.pdf"), width = 16, height = 10)
            pheatmap(
                heatmap_mat,
                annotation_col = annotation_col,
                annotation_colors = list(Group = group_colors),
                color = viridis(100),
                cluster_rows = TRUE,
                cluster_cols = TRUE,
                fontsize_row = 6,
                fontsize_col = 8,
                main = sprintf(
                    "Top %d Significant KEGG Orthologs (p < 0.05)",
                    nrow(top30_sig)
                )
            )
            dev.off()
            cat(sprintf("  Heatmap saved (%d KOs).\n", nrow(top30_sig)))
        },
        error = function(e) {
            cat(sprintf("  WARNING: Heatmap generation failed (%s). Skipping.\n", e$message))
            try(dev.off(), silent = TRUE)
        }
    )
} else {
    cat("  No significant KOs for heatmap. Skipping.\n")
}

#-------------------------------------------------------------------------------
# BARPLOT OF SIGNIFICANT KOs  (with pairwise significance annotations)
#-------------------------------------------------------------------------------
if (nrow(sig_kos) > 0) {
    cat("Creating significant KO barplot...\n")

    # Take top 15 significant KOs (sorted by p-value)
    top_sig <- sig_kos %>% slice_head(n = 15)

    # Create short labels: KO ID + truncated description
    top_sig <- top_sig %>%
        mutate(KO_Label = if (exists("ko_desc")) {
            paste0(KO, ": ", substr(ko_desc[KO], 1, 35))
        } else {
            KO
        })

    # Pivot to long format for plotting
    plot_data <- top_sig %>%
        mutate(KO_Label = factor(KO_Label, levels = rev(KO_Label))) %>%
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

    # For each KO, run pairwise Wilcoxon on the raw relative abundances
    annot_list <- list()
    for (i in seq_len(nrow(top_sig))) {
        ko_id <- top_sig$KO[i]
        ko_label <- top_sig$KO_Label[i]
        ko_vals <- as.numeric(ko_rel[ko_id, ])
        max_abund <- max(
            top_sig$Mean_Control[i], top_sig$Mean_LC_COPD[i],
            top_sig$Mean_LC_woCOPD[i]
        )

        for (j in seq_along(comparisons)) {
            g1 <- comparisons[[j]][1]
            g2 <- comparisons[[j]][2]
            vals1 <- ko_vals[metadata$group == g1]
            vals2 <- ko_vals[metadata$group == g2]

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
                KO_Label = ko_label,
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
    annot_df$KO_Label <- factor(annot_df$KO_Label,
        levels = levels(plot_data$KO_Label)
    )

    # Position symbols using a GLOBAL offset so low-abundance KOs
    # don't have their symbols piled on top of each other.
    global_max <- max(annot_df$max_abund, na.rm = TRUE)
    step_size <- global_max * 0.08 # fixed gap between comparison columns

    annot_df <- annot_df %>%
        mutate(
            y_pos = max_abund + step_size * comp_idx
        )

    # Build the base barplot
    ko_bar <- ggplot(plot_data, aes(
        x = KO_Label,
        y = Abundance, fill = Group
    )) +
        geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
        coord_flip() +
        scale_fill_manual(values = group_colors) +
        labs(
            x = "", y = "Mean Relative Abundance (%)",
            title = "Top Differentially Abundant KOs (p < 0.05)"
        ) +
        theme_publication() +
        theme(axis.text.y = element_text(size = 7))

    # Add significance annotations
    ko_bar <- ko_bar +
        geom_text(
            data = annot_df,
            aes(
                x = KO_Label,
                y = y_pos,
                label = symbol
            ),
            inherit.aes = FALSE,
            size = 3.5, fontface = "bold", hjust = 0.5, vjust = 0.5
        )

    # Expand x-axis to make room for annotations
    ko_bar <- ko_bar +
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
    ko_bar <- ko_bar +
        labs(caption = sig_legend) +
        theme(plot.caption = element_text(hjust = 0, size = 7, face = "italic"))

    ggsave(file.path(output_dir, "significant_ko_barplot.pdf"),
        ko_bar,
        width = 14, height = 8, dpi = 300
    )
    ggsave(file.path(output_dir, "significant_ko_barplot.png"),
        ko_bar,
        width = 14, height = 8, dpi = 300
    )
    cat("  Barplot with pairwise annotations saved.\n")
}

#-------------------------------------------------------------------------------
# SAVE RESULTS
#-------------------------------------------------------------------------------
cat("\nSaving results...\n")

write.csv(kw_results, file.path(tables_dir, "ko_differential_analysis.csv"),
    row.names = FALSE
)

if (nrow(sig_kos) > 0) {
    write.csv(sig_kos, file.path(tables_dir, "ko_significant.csv"),
        row.names = FALSE
    )
}

# Save PERMANOVA result
permanova_df <- data.frame(
    Source = rownames(permanova),
    permanova
)
write.csv(permanova_df, file.path(tables_dir, "ko_permanova.csv"),
    row.names = FALSE
)

cat("\n============================================\n")
cat("KO Visualization Complete!\n")
cat("============================================\n")
