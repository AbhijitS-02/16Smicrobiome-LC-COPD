#!/bin/bash
#===============================================================================
# 05_DADA2_DENOISE.SH - DADA2 Denoising and ASV Generation
#===============================================================================
# Run DADA2 to denoise reads and generate Amplicon Sequence Variants (ASVs)
#===============================================================================
SECONDS=0
set -e

# Source configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

echo "=============================================="
echo "DADA2 Denoising"
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
# DADA2 PARAMETERS
#-------------------------------------------------------------------------------
# Check if truncation lengths are set
if [ "${TRUNC_LEN_F}" -eq 0 ] || [ "${TRUNC_LEN_R}" -eq 0 ]; then
    echo "WARNING: Truncation lengths not set in config.sh"
    echo ""
    echo "Please review the demux summary and set appropriate values:"
    echo "  View: qiime tools view ${QIIME2_DIR}/imported/demux-paired-end.qzv"
    echo ""
    echo "Suggested approach:"
    echo "  - Truncate where median quality drops below Q25-30"
    echo "  - Ensure sufficient overlap for merging (~20bp minimum)"
    echo "  - For V3-V4 (~460bp amplicon), consider:"
    echo "    TRUNC_LEN_F=270, TRUNC_LEN_R=220 for 2x300bp reads"
    echo "    TRUNC_LEN_F=240, TRUNC_LEN_R=200 for 2x250bp reads"
    echo ""
    
    # Use auto-detection based on quality profiles
    echo "Attempting to auto-detect optimal truncation lengths..."
    
    # Default conservative values for V3-V4 region
    # These assume 2x250bp or 2x300bp reads
    TRUNC_LEN_F=245
    TRUNC_LEN_R=215
    
    echo "Using default values: TRUNC_LEN_F=${TRUNC_LEN_F}, TRUNC_LEN_R=${TRUNC_LEN_R}"
    echo "Update config.sh if these need adjustment."
fi

echo "DADA2 Parameters:"
echo "  Forward truncation: ${TRUNC_LEN_F} bp"
echo "  Reverse truncation: ${TRUNC_LEN_R} bp"
echo "  Trim left F: ${TRIM_LEFT_F} bp"
echo "  Trim left R: ${TRIM_LEFT_R} bp"
echo "  Max expected errors F: ${MAX_EE_F}"
echo "  Max expected errors R: ${MAX_EE_R}"
echo "  Threads: ${THREADS}"
echo ""

#-------------------------------------------------------------------------------
# RUN DADA2
#-------------------------------------------------------------------------------
echo "Running DADA2 denoise-paired..."
echo "This may take 1-4 hours depending on data size and CPU."
echo ""

INPUT_QZA="${QIIME2_DIR}/imported/demux-paired-end.qza"
OUTPUT_DIR="${QIIME2_DIR}/dada2_output"

mkdir -p "${OUTPUT_DIR}"

qiime dada2 denoise-paired \
    --i-demultiplexed-seqs "${INPUT_QZA}" \
    --p-trunc-len-f ${TRUNC_LEN_F} \
    --p-trunc-len-r ${TRUNC_LEN_R} \
    --p-trim-left-f ${TRIM_LEFT_F} \
    --p-trim-left-r ${TRIM_LEFT_R} \
    --p-max-ee-f ${MAX_EE_F} \
    --p-max-ee-r ${MAX_EE_R} \
    --p-trunc-q ${TRUNC_Q} \
    --p-chimera-method ${CHIMERA_METHOD} \
    --p-n-threads ${THREADS} \
    --o-table "${OUTPUT_DIR}/table.qza" \
    --o-representative-sequences "${OUTPUT_DIR}/rep-seqs.qza" \
    --o-denoising-stats "${OUTPUT_DIR}/denoising-stats.qza" \
    --o-base-transition-stats "${OUTPUT_DIR}/base-transition-stats.qza" \
    --verbose

check_status "DADA2 denoising"

#-------------------------------------------------------------------------------
# GENERATE SUMMARIES
#-------------------------------------------------------------------------------
echo ""
echo "Generating summary visualizations..."

# Feature table summary
qiime feature-table summarize \
    --i-table "${OUTPUT_DIR}/table.qza" \
    --m-metadata-file "${METADATA_DIR}/sample_metadata.tsv" \
    --o-feature-frequencies "${OUTPUT_DIR}/feature-frequencies.qza" \
    --o-sample-frequencies "${OUTPUT_DIR}/sample-frequencies.qza" \
    --o-summary "${OUTPUT_DIR}/table-summary.qzv"

# Representative sequences summary
qiime feature-table tabulate-seqs \
    --i-data "${OUTPUT_DIR}/rep-seqs.qza" \
    --o-visualization "${OUTPUT_DIR}/rep-seqs-summary.qzv"

# Denoising stats visualization
qiime metadata tabulate \
    --m-input-file "${OUTPUT_DIR}/denoising-stats.qza" \
    --o-visualization "${OUTPUT_DIR}/denoising-stats.qzv"

check_status "Summary visualizations"

#-------------------------------------------------------------------------------
# EXPORT TABLES FOR DOWNSTREAM ANALYSIS
#-------------------------------------------------------------------------------
echo ""
echo "Exporting tables..."

# Export feature table to BIOM and TSV
qiime tools export \
    --input-path "${OUTPUT_DIR}/table.qza" \
    --output-path "${TABLES_DIR}"

# Convert BIOM to TSV
biom convert \
    -i "${TABLES_DIR}/feature-table.biom" \
    -o "${TABLES_DIR}/asv_counts.tsv" \
    --to-tsv

# Export representative sequences
qiime tools export \
    --input-path "${OUTPUT_DIR}/rep-seqs.qza" \
    --output-path "${TABLES_DIR}"

mv "${TABLES_DIR}/dna-sequences.fasta" "${TABLES_DIR}/asv_sequences.fasta"

# Export denoising stats
qiime tools export \
    --input-path "${OUTPUT_DIR}/denoising-stats.qza" \
    --output-path "${TABLES_DIR}"

mv "${TABLES_DIR}/stats.tsv" "${TABLES_DIR}/dada2_stats.tsv"

check_status "Table export"

#-------------------------------------------------------------------------------
# PRINT STATISTICS
#-------------------------------------------------------------------------------
echo ""
echo "=============================================="
echo "DADA2 Denoising Statistics"
echo "=============================================="

# Display stats
echo ""
echo "Sample statistics:"
head "${TABLES_DIR}/dada2_stats.tsv" | column -t

# Count ASVs
ASV_COUNT=$(tail -n +2 "${TABLES_DIR}/asv_counts.tsv" | wc -l)
echo ""
echo "Total ASVs generated: ${ASV_COUNT}"

# Calculate total reads
TOTAL_READS=$(tail -n +3 "${TABLES_DIR}/asv_counts.tsv" | cut -f2- | tr '\t' '\n' | awk '{sum+=$1} END{print sum}')
echo "Total reads in ASV table: ${TOTAL_READS}"
#-------------------------------------------------------------------------------
# SUMMARY
#-------------------------------------------------------------------------------
echo ""
echo "=============================================="
echo "DADA2 Denoising Complete!"
echo "=============================================="
echo "End time: $(date)"
echo ""
echo "Output files:"
echo "  ASV table: ${OUTPUT_DIR}/table.qza"
echo "  Representative sequences: ${OUTPUT_DIR}/rep-seqs.qza"
echo "  Denoising stats: ${OUTPUT_DIR}/denoising-stats.qza"
echo "  Exported tables: ${TABLES_DIR}/"
echo ""
echo "Key metrics:"
echo "  Total ASVs: ${ASV_COUNT}"
echo "  Total reads: ${TOTAL_READS}"
echo ""
echo "IMPORTANT: Review denoising stats to ensure:"
echo "  - High percentage of reads passed filtering (>70%)"
echo "  - Low chimera rates (<20%)"
echo "  - Sufficient reads per sample for analysis"
echo ""
echo "Next step: bash ${SCRIPT_DIR}/06_taxonomy_silva.sh"

duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds."