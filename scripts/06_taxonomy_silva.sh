#!/bin/bash
#===============================================================================
# 06_TAXONOMY_SILVA.SH - Taxonomic Classification with SILVA Database
#===============================================================================
# Classify ASVs using the SILVA 138.1 database
#===============================================================================
SECONDS=0
set -e

# Source configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

echo "=============================================="
echo "Taxonomic Classification - SILVA ${SILVA_VERSION}"
echo "=============================================="
echo "Start time: $(date)"
echo ""

# Check if QIIME2 environment is active
if [[ "${CONDA_DEFAULT_ENV}" != "${QIIME2_ENV}" ]]; then
    echo "ERROR: QIIME2 environment is not active."
    echo "Expected environment: ${QIIME2_ENV}"
    echo "Current environment: ${CONDA_DEFAULT_ENV:-None}"
    echo ""
    echo "Please activate the environment before running this script:"
    echo "  conda activate ${QIIME2_ENV}"
    exit 1
fi

echo "Environment confirmed: ${CONDA_DEFAULT_ENV}"
echo ""

#-------------------------------------------------------------------------------
# CHECK CLASSIFIER
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# PREPARE CLASSIFIER
#-------------------------------------------------------------------------------
# Define training directory
SILVA_TRAINING_DIR="${DATABASE_DIR}/silva/training"
mkdir -p "${SILVA_TRAINING_DIR}"

if [ ! -f "${SILVA_CLASSIFIER}" ]; then
    echo "Classifier not found at: ${SILVA_CLASSIFIER}"
    echo "Initiating download and training process..."
    echo ""

    # 1. Reference sequences
    SILVA_SEQS="${SILVA_TRAINING_DIR}/silva-138-99-seqs.qza"
    if [ ! -f "${SILVA_SEQS}" ]; then
        echo "Downloading Reference sequences..."
        wget -q --show-progress -O "${SILVA_SEQS}" "${SILVA_SEQS_URL}"
    else
        echo "Reference sequences found: ${SILVA_SEQS}"
    fi

    # 2. Reference taxonomy
    SILVA_TAX="${SILVA_TRAINING_DIR}/silva-138-99-tax.qza"
    if [ ! -f "${SILVA_TAX}" ]; then
        echo "Downloading Reference taxonomy..."
        wget -q --show-progress -O "${SILVA_TAX}" "${SILVA_TAX_URL}"
    else
        echo "Reference taxonomy found: ${SILVA_TAX}"
    fi

    # 3. Extract reads
    # Using user-specified primers (515F/806R)
    SILVA_EXTRACTED="${SILVA_TRAINING_DIR}/silva-138-99-515-806-seqs.qza"
    if [ ! -f "${SILVA_EXTRACTED}" ]; then
        echo "Extracting reads for classifier training..."
        qiime feature-classifier extract-reads \
          --i-sequences "${SILVA_SEQS}" \
          --p-f-primer GTGYCAGCMGCCGCGGTAA \
          --p-r-primer GGACTACNVGGGTWTCTAAT \
          --p-trunc-len 0 \
          --o-reads "${SILVA_EXTRACTED}"
    else
        echo "Extracted reads found: ${SILVA_EXTRACTED}"
    fi

    # Refitting of classifier is required as a lower version of scikit-learn was used to create the 
    # qiime2 silva training artifact, whereas the recent qiime2 version has the upgraded version of
    # scikit-learn which is not compatible with the older artifact.

    # Around 40 GB of RAM (12 GB alloted from system and 28 GB from swap) was required for the below step.
    # On a Ryzen 5 8645HS system with 6 cores and 12 threads, the below step took ~ 3 hours to run on WSL2.
    
    # 4. Fit classifier
    echo "Training Naive Bayes classifier (this may take a while)..."
    qiime feature-classifier fit-classifier-naive-bayes \
      --i-reference-reads "${SILVA_EXTRACTED}" \
      --i-reference-taxonomy "${SILVA_TAX}" \
      --o-classifier "${SILVA_CLASSIFIER}"
      
    check_status "SILVA Classifier training"
else
    echo "Found existing classifier: ${SILVA_CLASSIFIER}"
fi

echo "Using classifier: ${SILVA_CLASSIFIER}"
echo ""

#-------------------------------------------------------------------------------
# CLASSIFY SEQUENCES
#-------------------------------------------------------------------------------
echo "Classifying ASV sequences with SILVA..."

INPUT_SEQS="${QIIME2_DIR}/dada2_output/rep-seqs.qza"
OUTPUT_DIR="${QIIME2_DIR}/taxonomy_silva"

mkdir -p "${OUTPUT_DIR}"

qiime feature-classifier classify-sklearn \
    --i-classifier "${SILVA_CLASSIFIER}" \
    --i-reads "${INPUT_SEQS}" \
    --p-n-jobs 1 \
    --p-confidence 0.7 \
    --p-read-orientation 'same' \
    --o-classification "${OUTPUT_DIR}/taxonomy.qza"

check_status "SILVA classification"

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
# FILTER OUT NON-BACTERIAL/ARCHAEAL SEQUENCES
#-------------------------------------------------------------------------------
echo ""
echo "Filtering taxonomy..."

# Remove mitochondria and chloroplast
qiime taxa filter-table \
    --i-table "${QIIME2_DIR}/dada2_output/table.qza" \
    --i-taxonomy "${OUTPUT_DIR}/taxonomy.qza" \
    --p-exclude "mitochondria,chloroplast" \
    --o-filtered-table "${OUTPUT_DIR}/table-filtered.qza"

# Remove unassigned at kingdom level
qiime taxa filter-table \
    --i-table "${OUTPUT_DIR}/table-filtered.qza" \
    --i-taxonomy "${OUTPUT_DIR}/taxonomy.qza" \
    --p-include "d__Bacteria,d__Archaea" \
    --o-filtered-table "${OUTPUT_DIR}/table-bacteria-archaea.qza"

# Summary of filtered table
qiime feature-table summarize \
  --i-table "${OUTPUT_DIR}/table-bacteria-archaea.qza" \
  --m-metadata-file "${METADATA_DIR}/sample_metadata.tsv" \
  --o-feature-frequencies "${OUTPUT_DIR}/feature-frequencies.qza" \
  --o-sample-frequencies "${OUTPUT_DIR}/sample-frequencies.qza" \
  --o-summary "${OUTPUT_DIR}/table-filtered-summary.qzv"

check_status "Taxonomy filtering"

#-------------------------------------------------------------------------------
# COLLAPSE TO TAXONOMIC LEVELS
#-------------------------------------------------------------------------------
echo ""
echo "Collapsing to taxonomic levels..."

LEVELS=("2:Phylum" "3:Class" "4:Order" "5:Family" "6:Genus" "7:Species")

for LEVEL_INFO in "${LEVELS[@]}"; do
    LEVEL=$(echo "${LEVEL_INFO}" | cut -d: -f1)
    LEVEL_NAME=$(echo "${LEVEL_INFO}" | cut -d: -f2)
    
    echo "  Collapsing to level ${LEVEL} (${LEVEL_NAME})..."
    
    qiime taxa collapse \
        --i-table "${OUTPUT_DIR}/table-bacteria-archaea.qza" \
        --i-taxonomy "${OUTPUT_DIR}/taxonomy.qza" \
        --p-level ${LEVEL} \
        --o-collapsed-table "${OUTPUT_DIR}/table-level-${LEVEL}-${LEVEL_NAME}.qza"
done

check_status "Taxonomy collapse"

#-------------------------------------------------------------------------------
# EXPORT DATA
#-------------------------------------------------------------------------------
echo ""
echo "Exporting taxonomy data..."

# Export taxonomy
qiime tools export \
    --input-path "${OUTPUT_DIR}/taxonomy.qza" \
    --output-path "${TABLES_DIR}/silva"

mv "${TABLES_DIR}/silva/taxonomy.tsv" "${TABLES_DIR}/taxonomy_silva.tsv"

# Export filtered table
qiime tools export \
    --input-path "${OUTPUT_DIR}/table-bacteria-archaea.qza" \
    --output-path "${TABLES_DIR}/silva"

biom convert \
    -i "${TABLES_DIR}/silva/feature-table.biom" \
    -o "${TABLES_DIR}/asv_counts_filtered_silva.tsv" \
    --to-tsv

# Clean up
rm -rf "${TABLES_DIR}/silva"

check_status "Data export"

#-------------------------------------------------------------------------------
# TAXONOMY SUMMARY
#-------------------------------------------------------------------------------
echo ""
echo "=============================================="
echo "SILVA Taxonomy Summary"
echo "=============================================="

# Count classifications at each level
declare -A LEVEL_MAP=(
  [Kingdom]="d__"
  [Phylum]="p__"
  [Class]="c__"
  [Order]="o__"
  [Family]="f__"
  [Genus]="g__"
  [Species]="s__"
)

echo ""
echo "Classification completeness:"
for LEVEL in "${!LEVEL_MAP[@]}"; do
  PREFIX=${LEVEL_MAP[$LEVEL]}
  COUNT=$(grep -v "^Feature" "${TABLES_DIR}/taxonomy_silva.tsv" | \
          cut -f2 | grep -c "${PREFIX}[^;]*[a-zA-Z]")
  echo "  ${LEVEL}: ${COUNT} ASVs"
done


# Top phyla
echo ""
echo "Top 10 Phyla:"
cut -f2 "${TABLES_DIR}/taxonomy_silva.tsv" | \
    grep -v "Taxon" | \
    cut -d';' -f1-2 | \
    sort | uniq -c | sort -rn | head -10

#-------------------------------------------------------------------------------
# SUMMARY
#-------------------------------------------------------------------------------
echo ""
echo "=============================================="
echo "SILVA Classification Complete!"
echo "=============================================="
echo "End time: $(date)"
echo ""
echo "Output files:"
echo "  Taxonomy: ${OUTPUT_DIR}/taxonomy.qza"
echo "  Barplot: ${OUTPUT_DIR}/taxa-barplot.qzv"
echo "  Filtered table: ${OUTPUT_DIR}/table-bacteria-archaea.qza"
echo "  Collapsed tables: ${OUTPUT_DIR}/table-level-*.qza"
echo "  Exported: ${TABLES_DIR}/taxonomy_silva.tsv"
echo ""
echo "Next step: bash ${SCRIPT_DIR}/07_taxonomy_gtdb.sh"

duration=$SECONDS
echo ""
echo "Total time: $(($duration / 60)) minutes and $(($duration % 60)) seconds"