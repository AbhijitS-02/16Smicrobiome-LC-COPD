#!/bin/bash
#===============================================================================
# 03_FASTQC_POSTTRIM.SH - Post-trimming Quality Control
#===============================================================================
# Run FastQC on trimmed reads to verify quality improvement
#===============================================================================
SECONDS=0
set -e

# Source configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

echo "=============================================="
echo "Post-trimming Quality Control (FastQC)"
echo "=============================================="
echo "Start time: $(date)"
echo ""

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate "${QIIME2_ENV}"

#-------------------------------------------------------------------------------
# RUN FASTQC ON TRIMMED READS
#-------------------------------------------------------------------------------
OUTPUT_DIR="${QC_DIR}/posttrim_fastqc"
mkdir -p "${OUTPUT_DIR}"

echo "Running FastQC on trimmed reads..."
echo "Input directory: ${TRIMMED_DIR}"
echo "Output directory: ${OUTPUT_DIR}"
echo "Using ${THREADS} threads"
echo ""

# Check if trimmed files exist
TRIMMED_COUNT=$(ls -1 ${TRIMMED_DIR}/*_trimmed.fastq.gz 2>/dev/null | wc -l)
if [ "${TRIMMED_COUNT}" -eq 0 ]; then
    echo "ERROR: No trimmed FASTQ files found in ${TRIMMED_DIR}"
    echo "Please run 02_trim_primers.sh first."
    exit 1
fi

echo "Processing ${TRIMMED_COUNT} trimmed FASTQ files..."

# Run FastQC
fastqc \
    --outdir "${OUTPUT_DIR}" \
    --threads "${THREADS}" \
    --quiet \
    ${TRIMMED_DIR}/*_trimmed.fastq.gz

check_status "FastQC analysis"

#-------------------------------------------------------------------------------
# RUN MULTIQC
#-------------------------------------------------------------------------------
echo ""
echo "Aggregating results with MultiQC..."

MULTIQC_DIR="${QC_DIR}/multiqc/posttrim"
mkdir -p "${MULTIQC_DIR}"

multiqc \
    --outdir "${MULTIQC_DIR}" \
    --filename "posttrim_multiqc_report" \
    --title "Post-trimming QC Report - LC-COPD Microbiome Study" \
    --force \
    "${OUTPUT_DIR}"

check_status "MultiQC aggregation"

#-------------------------------------------------------------------------------
# COMPARE PRE VS POST TRIMMING
#-------------------------------------------------------------------------------
echo ""
echo "Generating comparison report..."

COMPARISON_DIR="${QC_DIR}/multiqc/comparison"
mkdir -p "${COMPARISON_DIR}"

# Run MultiQC on both pre and post trimming
multiqc \
    --outdir "${COMPARISON_DIR}" \
    --filename "pretrim_vs_posttrim_comparison" \
    --title "QC Comparison: Pre vs Post Trimming" \
    --force \
    "${QC_DIR}/pretrim_fastqc" \
    "${QC_DIR}/posttrim_fastqc"

check_status "Comparison report"

#-------------------------------------------------------------------------------
# QUALITY SUMMARY
#-------------------------------------------------------------------------------
echo ""
echo "Creating quality summary..."

SUMMARY_FILE="${OUTPUT_DIR}/posttrim_quality_summary.txt"

echo "Post-trimming Quality Summary" > "${SUMMARY_FILE}"
echo "Generated: $(date)" >> "${SUMMARY_FILE}"
echo "==============================" >> "${SUMMARY_FILE}"
echo "" >> "${SUMMARY_FILE}"

# Count reads in trimmed files
echo "Read counts after trimming:" >> "${SUMMARY_FILE}"
echo "" >> "${SUMMARY_FILE}"

TOTAL_READS=0
for file in ${TRIMMED_DIR}/*_R1_trimmed.fastq.gz; do
    SAMPLE=$(basename "${file}" _R1_trimmed.fastq.gz)
    READS=$(zcat "${file}" | wc -l)
    READS=$((READS / 4))  # FASTQ has 4 lines per read
    TOTAL_READS=$((TOTAL_READS + READS))
    echo "  ${SAMPLE}: ${READS} read pairs" >> "${SUMMARY_FILE}"
done

echo "" >> "${SUMMARY_FILE}"
echo "Total reads across all samples: ${TOTAL_READS}" >> "${SUMMARY_FILE}"
echo "Average reads per sample: $((TOTAL_READS / 30))" >> "${SUMMARY_FILE}"

conda deactivate

#-------------------------------------------------------------------------------
# SUMMARY
#-------------------------------------------------------------------------------
echo ""
echo "=============================================="
echo "Post-trimming QC Complete!"
echo "=============================================="
echo "End time: $(date)"
echo ""
echo "Output files:"
echo "  FastQC reports: ${OUTPUT_DIR}/"
echo "  MultiQC report: ${MULTIQC_DIR}/posttrim_multiqc_report.html"
echo "  Comparison report: ${COMPARISON_DIR}/pretrim_vs_posttrim_comparison.html"
echo "  Summary: ${SUMMARY_FILE}"
echo ""
echo "REVIEW the comparison report to verify:"
echo "  1. Quality scores improved"
echo "  2. Adapter content removed"
echo "  3. Read lengths appropriate for merging"
echo ""
echo "Next step: bash ${SCRIPT_DIR}/04_qiime2_import.sh"

duration=$SECONDS
echo "$((duration / 60)) minutes $((duration % 60)) seconds elapsed."