#!/bin/bash
#===============================================================================
# 09_DIVERSITY.SH - Alpha and Beta Diversity Analysis
#===============================================================================
# Calculate diversity metrics and perform statistical tests
#===============================================================================
SECONDS=0
set -e

# Source configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

echo "=============================================="
echo "Diversity Analysis"
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

# Use filtered table from SILVA classification
INPUT_TABLE="${QIIME2_DIR}/taxonomy_silva/table-bacteria-archaea.qza"
INPUT_TREE="${QIIME2_DIR}/phylogeny/rooted-tree.qza"
OUTPUT_DIR="${QIIME2_DIR}/diversity"

mkdir -p "${OUTPUT_DIR}"

#-------------------------------------------------------------------------------
# DETERMINE RAREFACTION DEPTH
#-------------------------------------------------------------------------------
echo "Determining appropriate rarefaction depth..."

# Get feature table summary
qiime feature-table summarize \
    --i-table "${INPUT_TABLE}" \
    --o-feature-frequencies "${OUTPUT_DIR}/feature-frequencies.qza" \
    --o-sample-frequencies "${OUTPUT_DIR}/sample-frequencies.qza" \
    --o-summary "${OUTPUT_DIR}/table-summary-for-rarefaction.qzv"


# Export to get min reads per sample
qiime tools export \
    --input-path "${OUTPUT_DIR}/sample-frequencies.qza" \
    --output-path "${OUTPUT_DIR}/table-summary-export"

# Parse minimum sample depth
if [ -f "${OUTPUT_DIR}/table-summary-export/metadata.tsv" ]; then
    MIN_DEPTH=$(awk -F'\t' 'NR>2 {gsub(",", "", $2); gsub(/\.0$/, "", $2); print $2}' \
            "${OUTPUT_DIR}/table-summary-export/metadata.tsv" | \
            sort -n | head -1)

    RAREFACTION_DEPTH=$(echo "scale=0; ${MIN_DEPTH} * 0.9 / 1" | bc)
else
    echo "ERROR: Could not find exported sample frequencies file"
    exit 1
fi

echo "Minimum sample depth: ${MIN_DEPTH}"
echo "Rarefaction depth: ${RAREFACTION_DEPTH}"
echo ""

# Update config
echo "RAREFACTION_DEPTH=${RAREFACTION_DEPTH}" >> "${OUTPUT_DIR}/rarefaction_depth.txt"

#-------------------------------------------------------------------------------
# RAREFACTION CURVES
#-------------------------------------------------------------------------------
echo "Generating rarefaction curves..."

qiime diversity alpha-rarefaction \
    --i-table "${INPUT_TABLE}" \
    --i-phylogeny "${INPUT_TREE}" \
    --p-max-depth ${RAREFACTION_DEPTH} \
    --p-steps 20 \
    --m-metadata-file "${METADATA_DIR}/sample_metadata.tsv" \
    --o-visualization "${OUTPUT_DIR}/alpha-rarefaction.qzv"

check_status "Rarefaction curves"

#-------------------------------------------------------------------------------
# CORE DIVERSITY METRICS
#-------------------------------------------------------------------------------
echo ""
echo "Calculating core diversity metrics..."

qiime diversity core-metrics-phylogenetic \
    --i-table "${INPUT_TABLE}" \
    --i-phylogeny "${INPUT_TREE}" \
    --p-sampling-depth ${RAREFACTION_DEPTH} \
    --m-metadata-file "${METADATA_DIR}/sample_metadata.tsv" \
    --p-n-jobs-or-threads ${THREADS} \
    --output-dir "${OUTPUT_DIR}/core-metrics"

check_status "Core diversity metrics"

#-------------------------------------------------------------------------------
# ADDITIONAL ALPHA DIVERSITY METRICS
#-------------------------------------------------------------------------------
echo ""
echo "Calculating additional alpha diversity metrics..."

# Chao1
qiime diversity alpha \
    --i-table "${INPUT_TABLE}" \
    --p-metric chao1 \
    --o-alpha-diversity "${OUTPUT_DIR}/core-metrics/chao1_vector.qza"

# Simpson
qiime diversity alpha \
    --i-table "${INPUT_TABLE}" \
    --p-metric simpson \
    --o-alpha-diversity "${OUTPUT_DIR}/core-metrics/simpson_vector.qza"

# Fisher's alpha (if available)
qiime diversity alpha \
    --i-table "${INPUT_TABLE}" \
    --p-metric fisher_alpha \
    --o-alpha-diversity "${OUTPUT_DIR}/core-metrics/fisher_alpha_vector.qza" 2>/dev/null || \
    echo "  Fisher's alpha not available"

# ACE
qiime diversity alpha \
    --i-table "${INPUT_TABLE}" \
    --p-metric ace \
    --o-alpha-diversity "${OUTPUT_DIR}/core-metrics/ace_vector.qza" 2>/dev/null || \
    echo "  ACE not available"

check_status "Additional alpha metrics"

#-------------------------------------------------------------------------------
# ALPHA DIVERSITY STATISTICAL TESTS
#-------------------------------------------------------------------------------
echo ""
echo "Running alpha diversity statistical tests..."

# List of alpha metrics to test
ALPHA_METRICS=("observed_features" "shannon" "pielou_e" "faith_pd" "chao1" "simpson")

for METRIC in "${ALPHA_METRICS[@]}"; do
    VECTOR="${OUTPUT_DIR}/core-metrics/${METRIC}_vector.qza"
    
    if [ -f "${VECTOR}" ]; then
        echo "  Testing ${METRIC}..."
        
        # Kruskal-Wallis (all groups)
        qiime diversity alpha-group-significance \
            --i-alpha-diversity "${VECTOR}" \
            --m-metadata-file "${METADATA_DIR}/sample_metadata.tsv" \
            --o-visualization "${OUTPUT_DIR}/${METRIC}-group-significance.qzv"
    fi
done

check_status "Alpha diversity statistics"

#-------------------------------------------------------------------------------
# BETA DIVERSITY STATISTICAL TESTS
#-------------------------------------------------------------------------------
echo ""
echo "Running beta diversity statistical tests..."

BETA_METRICS=("bray_curtis" "jaccard" "weighted_unifrac" "unweighted_unifrac")

for METRIC in "${BETA_METRICS[@]}"; do
    MATRIX="${OUTPUT_DIR}/core-metrics/${METRIC}_distance_matrix.qza"
    
    if [ -f "${MATRIX}" ]; then
        echo "  Testing ${METRIC}..."
        
        # PERMANOVA
        qiime diversity beta-group-significance \
            --i-distance-matrix "${MATRIX}" \
            --m-metadata-file "${METADATA_DIR}/sample_metadata.tsv" \
            --m-metadata-column group \
            --p-method permanova \
            --p-pairwise \
            --o-visualization "${OUTPUT_DIR}/${METRIC}-permanova.qzv"
        
        # ANOSIM
        qiime diversity beta-group-significance \
            --i-distance-matrix "${MATRIX}" \
            --m-metadata-file "${METADATA_DIR}/sample_metadata.tsv" \
            --m-metadata-column group \
            --p-method anosim \
            --p-pairwise \
            --o-visualization "${OUTPUT_DIR}/${METRIC}-anosim.qzv"
        
        # Beta dispersion test
        qiime diversity beta-group-significance \
            --i-distance-matrix "${MATRIX}" \
            --m-metadata-file "${METADATA_DIR}/sample_metadata.tsv" \
            --m-metadata-column group \
            --p-method permdisp \
            --p-pairwise \
            --o-visualization "${OUTPUT_DIR}/${METRIC}-permdisp.qzv"
    fi
done

check_status "Beta diversity statistics"

#-------------------------------------------------------------------------------
# PCoA EMPEROR PLOTS
#-------------------------------------------------------------------------------
echo ""
echo "Generating Emperor PCoA plots..."

for METRIC in "${BETA_METRICS[@]}"; do
    PCOA="${OUTPUT_DIR}/core-metrics/${METRIC}_pcoa_results.qza"
    
    if [ -f "${PCOA}" ]; then
        qiime emperor plot \
            --i-pcoa "${PCOA}" \
            --m-metadata-file "${METADATA_DIR}/sample_metadata.tsv" \
            --o-visualization "${OUTPUT_DIR}/${METRIC}-emperor.qzv"
    fi
done

check_status "Emperor plots"

#-------------------------------------------------------------------------------
# EXPORT DIVERSITY RESULTS
#-------------------------------------------------------------------------------
echo ""
echo "Exporting diversity results..."

mkdir -p "${TABLES_DIR}/diversity"

# Export alpha diversity vectors
for METRIC in "${ALPHA_METRICS[@]}"; do
    VECTOR="${OUTPUT_DIR}/core-metrics/${METRIC}_vector.qza"
    
    if [ -f "${VECTOR}" ]; then
        qiime tools export \
            --input-path "${VECTOR}" \
            --output-path "${TABLES_DIR}/diversity"
        
        mv "${TABLES_DIR}/diversity/alpha-diversity.tsv" \
           "${TABLES_DIR}/diversity/alpha_${METRIC}.tsv"
    fi
done

# Export beta diversity matrices
for METRIC in "${BETA_METRICS[@]}"; do
    MATRIX="${OUTPUT_DIR}/core-metrics/${METRIC}_distance_matrix.qza"
    
    if [ -f "${MATRIX}" ]; then
        qiime tools export \
            --input-path "${MATRIX}" \
            --output-path "${TABLES_DIR}/diversity"
        
        mv "${TABLES_DIR}/diversity/distance-matrix.tsv" \
           "${TABLES_DIR}/diversity/beta_${METRIC}.tsv"
    fi
done

check_status "Diversity export"

#-------------------------------------------------------------------------------
# SUMMARY
#-------------------------------------------------------------------------------
echo ""
echo "=============================================="
echo "Diversity Analysis Complete!"
echo "=============================================="
echo "End time: $(date)"
echo ""
echo "Rarefaction depth used: ${RAREFACTION_DEPTH}"
echo ""
echo "Output files:"
echo "  Core metrics: ${OUTPUT_DIR}/core-metrics/"
echo "  Rarefaction: ${OUTPUT_DIR}/alpha-rarefaction.qzv"
echo "  Statistical tests: ${OUTPUT_DIR}/*-significance.qzv"
echo "  Emperor plots: ${OUTPUT_DIR}/*-emperor.qzv"
echo "  Exported tables: ${TABLES_DIR}/diversity/"
echo ""
echo "Next step: bash ${SCRIPT_DIR}/10_picrust2.sh"

duration=$SECONDS
echo "Diversity analysis completed in $(($duration / 60)) minutes and $(($duration % 60)) seconds."
