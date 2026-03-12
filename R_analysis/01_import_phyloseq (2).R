#===============================================================================
# 01_IMPORT_PHYLOSEQ.R - Import QIIME2 Data into R/phyloseq
#===============================================================================
# Import ASV tables, taxonomy, tree, and metadata into phyloseq objects
#===============================================================================

# Set up environment
suppressPackageStartupMessages({
  library(phyloseq)
  library(qiime2R)
  library(tidyverse)
  library(yaml)
})

# Configuration
project_dir <- "D:/LC_COPD_microbiome"
qiime2_dir <- file.path(project_dir, "results/qiime2")
tables_dir <- file.path(project_dir, "results/tables")
metadata_file <- file.path(project_dir, "metadata/sample_metadata.tsv")
output_dir <- file.path(project_dir, "R_analysis/data")

# Create output directory
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("============================================\n")
cat("Importing QIIME2 Data into phyloseq\n")
cat("============================================\n\n")

#-------------------------------------------------------------------------------
# IMPORT METADATA
#-------------------------------------------------------------------------------
cat("Loading metadata...\n")

metadata <- read.table(metadata_file, 
                       header = TRUE, 
                       sep = "\t", 
                       row.names = 1,
                       stringsAsFactors = FALSE)

# Convert to factors with proper levels
metadata$group <- factor(metadata$group, 
                         levels = c("Control", "LC_COPD", "LC_woCOPD"))
metadata$disease_status <- factor(metadata$disease_status,
                                  levels = c("Healthy", "Lung_Cancer"))
metadata$copd_status <- factor(metadata$copd_status,
                               levels = c("No", "Yes"))

cat(sprintf("  Loaded %d samples\n", nrow(metadata)))
cat(sprintf("  Groups: %s\n", paste(levels(metadata$group), collapse = ", ")))

#-------------------------------------------------------------------------------
# FUNCTION TO CREATE PHYLOSEQ OBJECT
#-------------------------------------------------------------------------------
create_phyloseq_from_qiime2 <- function(table_file, taxonomy_file, tree_file, 
                                         metadata_df, name = "phyloseq") {
  
  cat(sprintf("\nCreating phyloseq object: %s\n", name))
  
  # Import ASV table
  cat("  Loading ASV table...\n")
  if (grepl("\\.qza$", table_file)) {
    asv_table <- read_qza(table_file)$data
  } else {
    # Read from TSV (skip first line which is a comment)
    asv_df <- read.table(table_file, header = TRUE, sep = "\t", 
                         skip = 1, row.names = 1, comment.char = "")
    asv_table <- as.matrix(asv_df)
  }
  
  # Import taxonomy
  cat("  Loading taxonomy...\n")
  if (grepl("\\.qza$", taxonomy_file)) {
    taxonomy_raw <- read_qza(taxonomy_file)$data
    # Parse taxonomy string
    tax_parsed <- taxonomy_raw %>%
      separate(Taxon, into = c("Kingdom", "Phylum", "Class", "Order", 
                               "Family", "Genus", "Species"),
               sep = "; ", fill = "right", remove = FALSE) %>%
      mutate(across(Kingdom:Species, ~gsub("^[a-z]__", "", .))) %>%
      mutate(across(Kingdom:Species, ~na_if(., "")))
    
    tax_matrix <- as.matrix(tax_parsed[, c("Kingdom", "Phylum", "Class", 
                                           "Order", "Family", "Genus", "Species")])
    rownames(tax_matrix) <- taxonomy_raw$Feature.ID
  } else {
    tax_df <- read.table(taxonomy_file, header = TRUE, sep = "\t", 
                         row.names = 1, stringsAsFactors = FALSE)
    # Parse taxonomy string
    tax_split <- strsplit(tax_df$Taxon, "; ")
    tax_matrix <- do.call(rbind, lapply(tax_split, function(x) {
      ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
      vals <- sapply(ranks, function(r) {
        idx <- grep(paste0("^", substr(r, 1, 1), "__"), x, ignore.case = TRUE)
        if (length(idx) > 0) gsub("^[a-z]__", "", x[idx[1]]) else NA
      })
      return(vals)
    }))
    rownames(tax_matrix) <- rownames(tax_df)
  }
  
  # Import phylogenetic tree
  cat("  Loading phylogenetic tree...\n")
  if (grepl("\\.qza$", tree_file)) {
    tree <- read_qza(tree_file)$data
  } else {
    tree <- ape::read.tree(tree_file)
  }
  
  # Create phyloseq object
  cat("  Assembling phyloseq object...\n")
  
  # Ensure sample names match
  common_samples <- intersect(colnames(asv_table), rownames(metadata_df))
  cat(sprintf("  Common samples: %d\n", length(common_samples)))
  
  asv_table <- asv_table[, common_samples]
  metadata_subset <- metadata_df[common_samples, , drop = FALSE]
  
  # Ensure ASV IDs match between table and taxonomy
  common_asvs <- intersect(rownames(asv_table), rownames(tax_matrix))
  cat(sprintf("  Common ASVs: %d\n", length(common_asvs)))
  
  asv_table <- asv_table[common_asvs, ]
  tax_matrix <- tax_matrix[common_asvs, ]
  
  # Build phyloseq
  ps <- phyloseq(
    otu_table(asv_table, taxa_are_rows = TRUE),
    tax_table(tax_matrix),
    sample_data(metadata_subset),
    phy_tree(tree)
  )
  
  cat(sprintf("  Final phyloseq: %d taxa, %d samples\n", ntaxa(ps), nsamples(ps)))
  
  return(ps)
}

#-------------------------------------------------------------------------------
# CREATE PHYLOSEQ WITH SILVA TAXONOMY
#-------------------------------------------------------------------------------
cat("\n--- SILVA Taxonomy ---\n")

ps_silva <- create_phyloseq_from_qiime2(
  table_file = file.path(qiime2_dir, "taxonomy_silva/table-bacteria-archaea.qza"),
  taxonomy_file = file.path(qiime2_dir, "taxonomy_silva/taxonomy.qza"),
  tree_file = file.path(qiime2_dir, "phylogeny/rooted-tree.qza"),
  metadata_df = metadata,
  name = "SILVA"
)

#-------------------------------------------------------------------------------
# CREATE PHYLOSEQ WITH GTDB TAXONOMY (if available)
#-------------------------------------------------------------------------------
gtdb_taxonomy <- file.path(qiime2_dir, "taxonomy_gtdb/taxonomy.qza")

if (file.exists(gtdb_taxonomy)) {
  cat("\n--- GTDB Taxonomy ---\n")
  
  ps_gtdb <- create_phyloseq_from_qiime2(
    table_file = file.path(qiime2_dir, "taxonomy_gtdb/table-filtered.qza"),
    taxonomy_file = gtdb_taxonomy,
    tree_file = file.path(qiime2_dir, "phylogeny/rooted-tree.qza"),
    metadata_df = metadata,
    name = "GTDB"
  )
} else {
  cat("\nGTDB taxonomy not found. Skipping.\n")
  ps_gtdb <- NULL
}

#-------------------------------------------------------------------------------
# BASIC FILTERING
#-------------------------------------------------------------------------------
cat("\n--- Filtering ---\n")

# Remove singletons (taxa appearing only once across all samples)
ps_silva_filt <- prune_taxa(taxa_sums(ps_silva) > 1, ps_silva)
cat(sprintf("After removing singletons: %d taxa\n", ntaxa(ps_silva_filt)))

# Remove taxa with less than 10 total count
ps_silva_filt <- prune_taxa(taxa_sums(ps_silva) > 10, ps_silva)
cat(sprintf("After removing taxa with less than 10 counts: %d taxa\n", ntaxa(ps_silva_filt)))

#-------------------------------------------------------------------------------
# SUMMARY STATISTICS
#-------------------------------------------------------------------------------
cat("\n============================================\n")
cat("Summary Statistics\n")
cat("============================================\n\n")

# Sample summary
sample_sums_df <- data.frame(
  sample = sample_names(ps_silva_filt),
  reads  = sample_sums(ps_silva_filt),
  group  = sample_data(ps_silva_filt)$group,
  disease_status = sample_data(ps_silva_filt)$disease_status,
  copd_status    = sample_data(ps_silva_filt)$copd_status,
  sample_type    = sample_data(ps_silva_filt)$sample_type
)

cat("Reads per group:\n")
sample_sums_df %>%
  group_by(group) %>%
  summarise(
    n = n(),
    mean_reads = round(mean(reads)),
    min_reads = min(reads),
    max_reads = max(reads),
    .groups = "drop"
  ) %>%
  print()

# Taxa summary by rank
cat("\nTaxa per taxonomic rank:\n")
tax_df <- as.data.frame(tax_table(ps_silva_filt))
for (rank in colnames(tax_df)) {
  n_unique <- length(unique(na.omit(tax_df[[rank]])))
  cat(sprintf("  %s: %d unique taxa\n", rank, n_unique))
}

#-------------------------------------------------------------------------------
# SAVE PHYLOSEQ OBJECTS
#-------------------------------------------------------------------------------
cat("\n--- Saving Objects ---\n")

# Save as RDS
saveRDS(ps_silva, file.path(output_dir, "phyloseq_silva_raw.rds"))
saveRDS(ps_silva_filt, file.path(output_dir, "phyloseq_silva_filtered.rds"))

if (!is.null(ps_gtdb)) {
  saveRDS(ps_gtdb, file.path(output_dir, "phyloseq_gtdb_raw.rds"))
}

# Save metadata
write.csv(sample_sums_df, file.path(output_dir, "sample_read_counts.csv"), 
          row.names = FALSE)

cat("\nFiles saved to:", output_dir, "\n")
cat("  - phyloseq_silva_raw.rds\n")
cat("  - phyloseq_silva_filtered.rds\n")
if (!is.null(ps_gtdb)) cat("  - phyloseq_gtdb_raw.rds\n")
cat("  - sample_read_counts.csv\n")

cat("\n============================================\n")
cat("Import Complete!\n")
cat("============================================\n")
