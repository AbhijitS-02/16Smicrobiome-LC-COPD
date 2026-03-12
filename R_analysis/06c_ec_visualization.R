# ===============================================================================
# 06c_EC_VISUALIZATION.R - PICRUSt2 EC (Enzyme Commission) Analysis Visualization
# ===============================================================================
# Mirrors the pathway visualization pipeline (06_picrust_visualization.R)
# but operates on Enzyme Commission (EC) metagenome predictions.
#
# Outputs:
#   - PCA of EC abundances (with PERMANOVA annotation)
#   - Heatmap of top 30 significant ECs (p < 0.05)
#   - Barplot of top 15 significant ECs with pairwise Wilcoxon annotations
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
cat("PICRUSt2 EC (Enzyme Commission) Analysis\n")
cat("============================================\n\n")

#-------------------------------------------------------------------------------
# LOAD DATA
#-------------------------------------------------------------------------------
cat("Loading EC data...\n")

# Load metadata
metadata <- read.table(metadata_file,
    header = TRUE, sep = "\t",
    row.names = 1, stringsAsFactors = FALSE
)
metadata$group <- factor(metadata$group, levels = c("Control", "LC_COPD", "LC_woCOPD"))

# Load EC metagenome predictions
ec_file <- file.path(tables_dir, "EC_metagenome_with_descriptions.tsv")

if (!file.exists(ec_file)) {
    ec_file <- file.path(picrust_dir, "EC_metagenome_with_descriptions.tsv")
}

if (!file.exists(ec_file)) {
    ec_file <- file.path(picrust_dir, "picrust2_out/EC_metagenome_out/pred_metagenome_unstrat.tsv")
    if (file.exists(paste0(ec_file, ".gz"))) {
        system(paste("gunzip -k", paste0(ec_file, ".gz")))
    }
}

if (file.exists(ec_file)) {
    ec_data <- read.table(ec_file,
        header = TRUE, sep = "\t",
        row.names = 1, check.names = FALSE,
        quote = "", comment.char = ""
    )
    cat(sprintf(
        "  Loaded %d ECs across %d columns\n",
        nrow(ec_data), ncol(ec_data)
    ))
} else {
    stop("EC file not found. Please run PICRUSt2 first.")
}

#-------------------------------------------------------------------------------
# PROCESS EC DATA
#-------------------------------------------------------------------------------
cat("\nProcessing EC data...\n")

# Separate description column if present
if ("description" %in% colnames(ec_data)) {
    ec_desc <- ec_data$description
    names(ec_desc) <- rownames(ec_data)
    ec_data <- ec_data[, !colnames(ec_data) %in% "description"]
}

# Ensure samples match metadata
common_samples <- intersect(colnames(ec_data), rownames(metadata))

cat(sprintf(
    "  EC columns (first 5): %s\n",
    paste(head(colnames(ec_data), 5), collapse = ", ")
))
cat(sprintf(
    "  Metadata rows (first 5):   %s\n",
    paste(head(rownames(metadata), 5), collapse = ", ")
))
cat(sprintf("  Matching samples: %d\n", length(common_samples)))

if (length(common_samples) == 0) {
    stop(
        "No matching sample IDs between EC data columns and metadata rows!\n",
        "  Check that sample IDs are consistent between files."
    )
}

ec_data <- ec_data[, common_samples]
metadata <- metadata[common_samples, , drop = FALSE]

# Convert to numeric (precaution)
ec_data <- as.data.frame(lapply(ec_data, as.numeric))
rownames(ec_data) <- names(ec_desc) # restore rownames after lapply

# Convert to relative abundance
ec_rel <- apply(ec_data, 2, function(x) x / sum(x) * 100)

cat(sprintf("  Final: %d ECs, %d samples\n", nrow(ec_rel), ncol(ec_rel)))

#-------------------------------------------------------------------------------
# EC PCA
#-------------------------------------------------------------------------------
cat("\nPerforming PCA on EC abundances...\n")

# Remove zero-variance rows before PCA (ECs with 0 in all samples)
ec_var_all <- apply(ec_rel, 1, var)
ec_rel_nzv <- ec_rel[ec_var_all > 0, ]
cat(sprintf("  Removed %d zero-variance ECs for PCA\n", sum(ec_var_all == 0)))

# Transpose for PCA
ec_t <- t(ec_rel_nzv)

# PCA
pca_result <- prcomp(ec_t, scale. = TRUE, center = TRUE)

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
        title = "PCA of Predicted EC Abundances",
        color = "Group"
    ) +
    theme_publication()

ggsave(file.path(output_dir, "ec_pca.pdf"),
    pca_plot,
    width = 8, height = 6, dpi = 300
)
ggsave(file.path(output_dir, "ec_pca.png"),
    pca_plot,
    width = 8, height = 6, dpi = 300
)

#-------------------------------------------------------------------------------
# PERMANOVA ON ECs
#-------------------------------------------------------------------------------
cat("Running PERMANOVA on EC abundances...\n")

# Calculate Bray-Curtis distance
ec_dist <- vegdist(ec_t, method = "bray")

# PERMANOVA
set.seed(42)
permanova <- adonis2(ec_dist ~ group, data = metadata, permutations = 200)

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

ggsave(file.path(output_dir, "ec_pca_annotated.pdf"),
    pca_plot_annotated,
    width = 8, height = 6, dpi = 300
)

#-------------------------------------------------------------------------------
# TOP DIFFERENTIAL ECs  (vectorized)
#-------------------------------------------------------------------------------
cat("\nIdentifying differentially abundant ECs...\n")

# Vectorized Kruskal-Wallis testing
p_values <- apply(ec_rel, 1, function(x) {
    kruskal.test(as.numeric(x) ~ metadata$group)$p.value
})

# Calculate means per group
group_idx <- split(seq_along(metadata$group), metadata$group)
means_control <- rowMeans(ec_rel[, group_idx[["Control"]], drop = FALSE])
means_lc_copd <- rowMeans(ec_rel[, group_idx[["LC_COPD"]], drop = FALSE])
means_lc_wocopd <- rowMeans(ec_rel[, group_idx[["LC_woCOPD"]], drop = FALSE])

# Build results data.frame
kw_results <- data.frame(
    EC = rownames(ec_rel),
    Description = if (exists("ec_desc")) {
        ec_desc[rownames(ec_rel)]
    } else {
        rownames(ec_rel)
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

# Significant ECs (using raw p-value)
sig_ecs <- kw_results %>% filter(p_value < 0.05)
cat(sprintf("  Significant ECs (p < 0.05): %d\n", nrow(sig_ecs)))

#-------------------------------------------------------------------------------
# HEATMAP OF TOP 30 SIGNIFICANT ECs  (p < 0.05)
#-------------------------------------------------------------------------------
cat("Creating EC heatmap...\n")

# Use top 30 significant ECs (already sorted by p-value)
if (nrow(sig_ecs) > 0) {
    top30_sig <- sig_ecs %>% slice_head(n = 30)
    top30_ecs <- ec_rel[top30_sig$EC, ]

    # Use description as row labels if available
    row_labels <- rownames(top30_ecs)
    if (exists("ec_desc")) {
        row_labels <- paste0(rownames(top30_ecs), ": ", ec_desc[rownames(top30_ecs)])
    }
    row_labels <- substr(row_labels, 1, 60) # Truncate long names

    # Annotation
    annotation_col <- data.frame(Group = metadata$group, row.names = rownames(metadata))

    # Log transform for visualization
    heatmap_mat <- log10(top30_ecs + 0.001)
    rownames(heatmap_mat) <- row_labels

    # Safe heatmap rendering via tryCatch
    tryCatch(
        {
            pdf(file.path(output_dir, "ec_heatmap_top30.pdf"), width = 16, height = 10)
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
                    "Top %d Significant Enzyme Commission Numbers (p < 0.05)",
                    nrow(top30_sig)
                )
            )
            dev.off()
            cat(sprintf("  Heatmap saved (%d ECs).\n", nrow(top30_sig)))
        },
        error = function(e) {
            cat(sprintf("  WARNING: Heatmap generation failed (%s). Skipping.\n", e$message))
            try(dev.off(), silent = TRUE)
        }
    )
} else {
    cat("  No significant ECs for heatmap. Skipping.\n")
}

#-------------------------------------------------------------------------------
# BARPLOT OF SIGNIFICANT ECs  (with pairwise significance annotations)
#-------------------------------------------------------------------------------
if (nrow(sig_ecs) > 0) {
    cat("Creating significant EC barplot...\n")

    # Take top 15 significant ECs (sorted by p-value)
    top_sig <- sig_ecs %>% slice_head(n = 15)

    # Create short labels: EC number + truncated description
    top_sig <- top_sig %>%
        mutate(EC_Label = if (exists("ec_desc")) {
            paste0(EC, ": ", substr(ec_desc[EC], 1, 35))
        } else {
            EC
        })

    # Pivot to long format for plotting
    plot_data <- top_sig %>%
        mutate(EC_Label = factor(EC_Label, levels = rev(EC_Label))) %>%
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

    # For each EC, run pairwise Wilcoxon on the raw relative abundances
    annot_list <- list()
    for (i in seq_len(nrow(top_sig))) {
        ec_id <- top_sig$EC[i]
        ec_label <- top_sig$EC_Label[i]
        ec_vals <- as.numeric(ec_rel[ec_id, ])
        max_abund <- max(
            top_sig$Mean_Control[i], top_sig$Mean_LC_COPD[i],
            top_sig$Mean_LC_woCOPD[i]
        )

        for (j in seq_along(comparisons)) {
            g1 <- comparisons[[j]][1]
            g2 <- comparisons[[j]][2]
            vals1 <- ec_vals[metadata$group == g1]
            vals2 <- ec_vals[metadata$group == g2]

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
                EC_Label = ec_label,
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
    annot_df$EC_Label <- factor(annot_df$EC_Label,
        levels = levels(plot_data$EC_Label)
    )

    # Position symbols using a GLOBAL offset so low-abundance ECs
    # don't have their symbols piled on top of each other.
    global_max <- max(annot_df$max_abund, na.rm = TRUE)
    step_size <- global_max * 0.08 # fixed gap between comparison columns

    annot_df <- annot_df %>%
        mutate(
            y_pos = max_abund + step_size * comp_idx
        )

    # Build the base barplot
    ec_bar <- ggplot(plot_data, aes(
        x = EC_Label,
        y = Abundance, fill = Group
    )) +
        geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
        coord_flip() +
        scale_fill_manual(values = group_colors) +
        labs(
            x = "", y = "Mean Relative Abundance (%)",
            title = "Top Differentially Abundant ECs (p < 0.05)"
        ) +
        theme_publication() +
        theme(axis.text.y = element_text(size = 7))

    # Add significance annotations
    ec_bar <- ec_bar +
        geom_text(
            data = annot_df,
            aes(
                x = EC_Label,
                y = y_pos,
                label = symbol
            ),
            inherit.aes = FALSE,
            size = 3.5, fontface = "bold", hjust = 0.5, vjust = 0.5
        )

    # Expand x-axis to make room for annotations
    ec_bar <- ec_bar +
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
    ec_bar <- ec_bar +
        labs(caption = sig_legend) +
        theme(plot.caption = element_text(hjust = 0, size = 7, face = "italic"))

    ggsave(file.path(output_dir, "significant_ec_barplot.pdf"),
        ec_bar,
        width = 14, height = 8, dpi = 300
    )
    ggsave(file.path(output_dir, "significant_ec_barplot.png"),
        ec_bar,
        width = 14, height = 8, dpi = 300
    )
    cat("  Barplot with pairwise annotations saved.\n")
}

#-------------------------------------------------------------------------------
# SAVE RESULTS
#-------------------------------------------------------------------------------
cat("\nSaving results...\n")

write.csv(kw_results, file.path(tables_dir, "ec_differential_analysis.csv"),
    row.names = FALSE
)

if (nrow(sig_ecs) > 0) {
    write.csv(sig_ecs, file.path(tables_dir, "ec_significant.csv"),
        row.names = FALSE
    )
}

# Save PERMANOVA result
permanova_df <- data.frame(
    Source = rownames(permanova),
    permanova
)
write.csv(permanova_df, file.path(tables_dir, "ec_permanova.csv"),
    row.names = FALSE
)

cat("\n============================================\n")
cat("EC Visualization Complete!\n")
cat("============================================\n")
