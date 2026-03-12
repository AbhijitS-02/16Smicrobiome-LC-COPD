#!/bin/bash
#===============================================================================
# 07_TAXONOMY_GTDB.SH - Taxonomic Classification with GTDB Database
#===============================================================================
# Classify ASVs using the GTDB r214 database
#===============================================================================

set -e

# Source configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

echo "=============================================="
echo "Taxonomic Classification - GTDB ${GTDB_VERSION}"
echo "=============================================="
echo "Start time: $(date)"
echo ""

# Activate QIIME2 environment
eval "$(conda shell.bash hook)"
conda activate "${QIIME2_ENV}"

#-------------------------------------------------------------------------------
# CHECK/CREATE GTDB CLASSIFIER
#-------------------------------------------------------------------------------
if [ ! -f "${GTDB_CLASSIFIER}" ]; then
    echo "GTDB classifier not found at: ${GTDB_CLASSIFIER}"
    echo ""
    echo "Option 1: Download pre-trained classifier (if available)"
    echo "Option 2: Train classifier using RESCRIPt (takes several hours)"
    echo ""
    echo "Attempting to use RESCRIPt to prepare GTDB classifier..."
    
    # Check if RESCRIPt is available
    if qiime rescript --help &>/dev/null; then
        echo "RESCRIPt plugin available. Training GTDB classifier..."
        
        mkdir -p "${DATABASE_DIR}/gtdb"
        
        # Download GTDB SSU sequences
        echo "Downloading GTDB SSU sequences..."
        qiime rescript get-gtdb-data \
            --p-version "214.1" \
            --p-domain "Bacteria" \
            --p-db-type "ssu-all" \
            --o-gtdb-taxonomy "${DATABASE_DIR}/gtdb/gtdb-taxonomy.qza" \
            --o-gtdb-sequences "${DATABASE_DIR}/gtdb/gtdb-sequences.qza" \
            --verbose || {
                echo "Warning: Could not download GTDB data via RESCRIPt"
                echo "Will use alternative VSEARCH-based classification"
            }
        
        if [ -f "${DATABASE_DIR}/gtdb/gtdb-sequences.qza" ]; then
            # Dereplicate sequences
            echo "Dereplicating sequences..."
            qiime rescript dereplicate \
                --i-sequences "${DATABASE_DIR}/gtdb/gtdb-sequences.qza" \
                --i-taxa "${DATABASE_DIR}/gtdb/gtdb-taxonomy.qza" \
                --p-mode 'uniq' \
                --o-dereplicated-sequences "${DATABASE_DIR}/gtdb/gtdb-seqs-derep.qza" \
                --o-dereplicated-taxa "${DATABASE_DIR}/gtdb/gtdb-tax-derep.qza"
            
            # Train classifier
            echo "Training classifier (this may take 1-2 hours)..."
            qiime feature-classifier fit-classifier-naive-bayes \
                --i-reference-reads "${DATABASE_DIR}/gtdb/gtdb-seqs-derep.qza" \
                --i-reference-taxonomy "${DATABASE_DIR}/gtdb/gtdb-tax-derep.qza" \
                --o-classifier "${GTDB_CLASSIFIER}"
            
            check_status "GTDB classifier training"
        fi
    else
        echo "RESCRIPt not available. Using VSEARCH-based classification instead."
    fi
fi

#-------------------------------------------------------------------------------
# CLASSIFY SEQUENCES
#-------------------------------------------------------------------------------
INPUT_SEQS="${QIIME2_DIR}/dada2_output/rep-seqs.qza"
OUTPUT_DIR="${QIIME2_DIR}/taxonomy_gtdb"

mkdir -p "${OUTPUT_DIR}"

if [ -f "${GTDB_CLASSIFIER}" ]; then
    echo "Classifying ASV sequences with GTDB..."
    
    qiime feature-classifier classify-sklearn \
        --i-classifier "${GTDB_CLASSIFIER}" \
        --i-reads "${INPUT_SEQS}" \
        --p-n-jobs ${THREADS} \
        --p-confidence 0.7 \
        --o-classification "${OUTPUT_DIR}/taxonomy.qza"
    
    check_status "GTDB classification"
    
elif [ -f "${DATABASE_DIR}/gtdb/gtdb-seqs-derep.qza" ]; then
    echo "Using VSEARCH consensus classification with GTDB..."
    
    qiime feature-classifier classify-consensus-vsearch \
        --i-query "${INPUT_SEQS}" \
        --i-reference-reads "${DATABASE_DIR}/gtdb/gtdb-seqs-derep.qza" \
        --i-reference-taxonomy "${DATABASE_DIR}/gtdb/gtdb-tax-derep.qza" \
        --p-perc-identity 0.97 \
        --p-threads ${THREADS} \
        --o-classification "${OUTPUT_DIR}/taxonomy.qza" \
        --o-search-results "${OUTPUT_DIR}/search-results.qza"
    
    check_status "GTDB VSEARCH classification"
else
    echo "ERROR: No GTDB classifier or reference database available."
    echo "Please either:"
    echo "  1. Download a pre-trained GTDB classifier"
    echo "  2. Install RESCRIPt plugin and re-run this script"
    echo ""
    echo "To install RESCRIPt:"
    echo "  conda install -c conda-forge -c bioconda -c qiime2 q2-rescript"
    echo ""
    echo "Skipping GTDB classification for now."
    
    conda deactivate
    exit 0
fi

#-------------------------------------------------------------------------------
# GENERATE VISUALIZATIONS
#-------------------------------------------------------------------------------
echo ""
echo "Generating taxonomy visualizations..."

# Taxonomy table
qiime metadata tabulate \
    --m-input-file "${OUTPUT_DIR}/taxonomy.qza" \
    --o-visualization "${OUTPUT_DIR}/taxonomy.qzv"

# Taxonomy barplot
qiime taxa barplot \
    --i-table "${QIIME2_DIR}/dada2_output/table.qza" \
    --i-taxonomy "${OUTPUT_DIR}/taxonomy.qza" \
    --m-metadata-file "${METADATA_DIR}/sample_metadata.tsv" \
    --o-visualization "${OUTPUT_DIR}/taxa-barplot.qzv"

check_status "Taxonomy visualizations"

#-------------------------------------------------------------------------------
# FILTER TAXONOMY
#-------------------------------------------------------------------------------
echo ""
echo "Filtering taxonomy..."

# Remove mitochondria and chloroplast
qiime taxa filter-table \
    --i-table "${QIIME2_DIR}/dada2_output/table.qza" \
    --i-taxonomy "${OUTPUT_DIR}/taxonomy.qza" \
    --p-exclude "mitochondria,chloroplast" \
    --o-filtered-table "${OUTPUT_DIR}/table-filtered.qza"

# Summary of filtered table
qiime feature-table summarize \
    --i-table "${OUTPUT_DIR}/table-filtered.qza" \
    --m-sample-metadata-file "${METADATA_DIR}/sample_metadata.tsv" \
    --o-visualization "${OUTPUT_DIR}/table-filtered-summary.qzv"

check_status "Taxonomy filtering"

#-------------------------------------------------------------------------------
# EXPORT DATA
#-------------------------------------------------------------------------------
echo ""
echo "Exporting taxonomy data..."

# Export taxonomy
qiime tools export \
    --input-path "${OUTPUT_DIR}/taxonomy.qza" \
    --output-path "${TABLES_DIR}/gtdb"

mv "${TABLES_DIR}/gtdb/taxonomy.tsv" "${TABLES_DIR}/taxonomy_gtdb.tsv"

# Export filtered table
qiime tools export \
    --input-path "${OUTPUT_DIR}/table-filtered.qza" \
    --output-path "${TABLES_DIR}/gtdb"

biom convert \
    -i "${TABLES_DIR}/gtdb/feature-table.biom" \
    -o "${TABLES_DIR}/asv_counts_filtered_gtdb.tsv" \
    --to-tsv

# Clean up
rm -rf "${TABLES_DIR}/gtdb"

check_status "Data export"

#-------------------------------------------------------------------------------
# COMPARE SILVA VS GTDB
#-------------------------------------------------------------------------------
echo ""
echo "=============================================="
echo "SILVA vs GTDB Classification Comparison"
echo "=============================================="

if [ -f "${TABLES_DIR}/taxonomy_silva.tsv" ] && [ -f "${TABLES_DIR}/taxonomy_gtdb.tsv" ]; then
    echo ""
    echo "SILVA classifications:"
    echo "  Total ASVs: $(tail -n +2 ${TABLES_DIR}/taxonomy_silva.tsv | wc -l)"
    
    echo ""
    echo "GTDB classifications:"
    echo "  Total ASVs: $(tail -n +2 ${TABLES_DIR}/taxonomy_gtdb.tsv | wc -l)"
    
    # Create comparison file
    echo "Creating comparison file..."
    paste \
        <(cut -f1,2 "${TABLES_DIR}/taxonomy_silva.tsv") \
        <(cut -f2 "${TABLES_DIR}/taxonomy_gtdb.tsv") \
        > "${TABLES_DIR}/taxonomy_comparison_silva_gtdb.tsv"
    
    # Add header
    sed -i '1s/.*/Feature_ID\tSILVA_Taxonomy\tGTDB_Taxonomy/' \
        "${TABLES_DIR}/taxonomy_comparison_silva_gtdb.tsv"
    
    echo "Comparison saved to: ${TABLES_DIR}/taxonomy_comparison_silva_gtdb.tsv"
fi

conda deactivate

#-------------------------------------------------------------------------------
# SUMMARY
#-------------------------------------------------------------------------------
echo ""
echo "=============================================="
echo "GTDB Classification Complete!"
echo "=============================================="
echo "End time: $(date)"
echo ""
echo "Output files:"
echo "  Taxonomy: ${OUTPUT_DIR}/taxonomy.qza"
echo "  Barplot: ${OUTPUT_DIR}/taxa-barplot.qzv"
echo "  Filtered table: ${OUTPUT_DIR}/table-filtered.qza"
echo "  Exported: ${TABLES_DIR}/taxonomy_gtdb.tsv"
echo ""
echo "Next step: bash ${SCRIPT_DIR}/08_phylogeny.sh"
