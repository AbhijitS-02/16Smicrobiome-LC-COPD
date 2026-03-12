#!/bin/bash
#===============================================================================
# 01_FASTQC_PRETRIM.SH - Pre-trimming Quality Control
#===============================================================================
# Run FastQC on raw reads to assess quality before trimming
#===============================================================================
SECONDS=0
set -e

# Source configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

echo "=============================================="
echo "Pre-trimming Quality Control (FastQC)"
echo "=============================================="
echo "Start time: $(date)"
echo ""

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate "${QIIME2_ENV}"

#-------------------------------------------------------------------------------
# RUN FASTQC
#-------------------------------------------------------------------------------
OUTPUT_DIR="${QC_DIR}/pretrim_fastqc"
mkdir -p "${OUTPUT_DIR}"

echo "Running FastQC on raw reads..."
echo "Input directory: ${RAW_DATA_DIR}"
echo "Output directory: ${OUTPUT_DIR}"
echo "Using ${THREADS} threads"
echo ""

# Count total files
TOTAL_FILES=$(ls -1 ${RAW_DATA_DIR}/*.fastq.gz 2>/dev/null | wc -l)
echo "Processing ${TOTAL_FILES} FASTQ files..."

# Run FastQC
fastqc \
    --outdir "${OUTPUT_DIR}" \
    --threads "${THREADS}" \
    --quiet \
    ${RAW_DATA_DIR}/*.fastq.gz

check_status "FastQC analysis"

#-------------------------------------------------------------------------------
# RUN MULTIQC
#-------------------------------------------------------------------------------
echo ""
echo "Aggregating results with MultiQC..."

MULTIQC_DIR="${QC_DIR}/multiqc/pretrim"
mkdir -p "${MULTIQC_DIR}"

multiqc \
    --outdir "${MULTIQC_DIR}" \
    --filename "pretrim_multiqc_report" \
    --title "Pre-trimming QC Report - LC-COPD Microbiome Study" \
    --force \
    "${OUTPUT_DIR}"

check_status "MultiQC aggregation"

#-------------------------------------------------------------------------------
# EXTRACT KEY METRICS
#-------------------------------------------------------------------------------
echo ""
echo "Extracting key quality metrics..."

# Create summary file
SUMMARY_FILE="${OUTPUT_DIR}/quality_summary.txt"

echo "Pre-trimming Quality Summary" > "${SUMMARY_FILE}"
echo "Generated: $(date)" >> "${SUMMARY_FILE}"
echo "==============================" >> "${SUMMARY_FILE}"
echo "" >> "${SUMMARY_FILE}"

# Parse FastQC data files for summary
for zip_file in ${OUTPUT_DIR}/*_fastqc.zip; do
    sample_name=$(basename "${zip_file}" _fastqc.zip)
    
    # Extract summary
    unzip -p "${zip_file}" "*/summary.txt" 2>/dev/null | while read status module; do
        if [ "${status}" == "FAIL" ]; then
            echo "${sample_name}: FAIL - ${module}" >> "${SUMMARY_FILE}"
        fi
    done
done

echo "" >> "${SUMMARY_FILE}"
echo "See MultiQC report for detailed analysis:" >> "${SUMMARY_FILE}"
echo "${MULTIQC_DIR}/pretrim_multiqc_report.html" >> "${SUMMARY_FILE}"

#-------------------------------------------------------------------------------
# DETERMINE READ LENGTHS
#-------------------------------------------------------------------------------
echo ""
echo "Determining read lengths from raw data..."

# Sample one file to get read length
SAMPLE_FILE=$(ls ${RAW_DATA_DIR}/*_R1.fastq.gz | head -1)
READ_LENGTH=$(zcat "${SAMPLE_FILE}" | head -2 | tail -1 | wc -c)
READ_LENGTH=$((READ_LENGTH - 1))  # Remove newline

echo "Detected read length: ${READ_LENGTH} bp"
echo "" >> "${SUMMARY_FILE}"
echo "Detected read length: ${READ_LENGTH} bp" >> "${SUMMARY_FILE}"

# Save read length for later use
echo "READ_LENGTH=${READ_LENGTH}" > "${OUTPUT_DIR}/read_length.txt"

conda deactivate

#-------------------------------------------------------------------------------
# SUMMARY
#-------------------------------------------------------------------------------
echo ""
echo "=============================================="
echo "Pre-trimming QC Complete!"
echo "=============================================="
echo "End time: $(date)"
echo ""
echo "Output files:"
echo "  FastQC reports: ${OUTPUT_DIR}/"
echo "  MultiQC report: ${MULTIQC_DIR}/pretrim_multiqc_report.html"
echo "  Summary: ${SUMMARY_FILE}"
echo ""
echo "IMPORTANT: Review the MultiQC report to:"
echo "  1. Confirm read lengths"
echo "  2. Check quality score distribution"
echo "  3. Identify adapter contamination"
echo "  4. Determine optimal trimming parameters"
echo ""
echo "Next step: bash ${SCRIPT_DIR}/02_trim_primers.sh"

duration=$SECONDS
echo "$((duration / 60)) minutes $((duration % 60)) seconds elapsed."
