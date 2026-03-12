#!/bin/bash
#===============================================================================
# 04_QIIME2_IMPORT.SH - Import Data into QIIME2
#===============================================================================
# Import trimmed paired-end reads into QIIME2 format
#===============================================================================
SECONDS=0
set -e

# Source configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

echo "=============================================="
echo "QIIME2 Data Import"
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

# Verify QIIME2 is working
qiime --version
echo ""

#-------------------------------------------------------------------------------
# CREATE MANIFEST FILE
#-------------------------------------------------------------------------------
echo "Creating manifest file..."

MANIFEST_FILE="${METADATA_DIR}/manifest.tsv"

# Create manifest header
echo -e "sample-id\tforward-absolute-filepath\treverse-absolute-filepath" > "${MANIFEST_FILE}"

# Add each sample
for R1 in ${TRIMMED_DIR}/*_R1_trimmed.fastq.gz; do
    SAMPLE=$(basename "${R1}" _R1_trimmed.fastq.gz)
    R2="${TRIMMED_DIR}/${SAMPLE}_R2_trimmed.fastq.gz"
    
    if [ -f "${R2}" ]; then
        echo -e "${SAMPLE}\t${R1}\t${R2}" >> "${MANIFEST_FILE}"
    else
        echo "WARNING: Missing R2 file for ${SAMPLE}"
    fi
done

echo "Manifest created: ${MANIFEST_FILE}"
echo "Samples in manifest: $(tail -n +2 ${MANIFEST_FILE} | wc -l)"

#-------------------------------------------------------------------------------
# IMPORT INTO QIIME2
#-------------------------------------------------------------------------------
echo ""
echo "Importing sequences into QIIME2..."

OUTPUT_QZA="${QIIME2_DIR}/imported/demux-paired-end.qza"

qiime tools import \
    --type 'SampleData[PairedEndSequencesWithQuality]' \
    --input-path "${MANIFEST_FILE}" \
    --output-path "${OUTPUT_QZA}" \
    --input-format PairedEndFastqManifestPhred33V2

check_status "QIIME2 import"

#-------------------------------------------------------------------------------
# GENERATE DEMUX SUMMARY
#-------------------------------------------------------------------------------
echo ""
echo "Generating demultiplex summary visualization..."

DEMUX_VIZ="${QIIME2_DIR}/imported/demux-paired-end.qzv"

qiime demux summarize \
    --i-data "${OUTPUT_QZA}" \
    --o-visualization "${DEMUX_VIZ}"

check_status "Demux summary"

# Export visualization for easy viewing
echo ""
echo "Exporting visualization data..."

DEMUX_EXPORT="${QIIME2_DIR}/imported/demux-summary-export"
qiime tools export \
    --input-path "${DEMUX_VIZ}" \
    --output-path "${DEMUX_EXPORT}"

#-------------------------------------------------------------------------------
# VALIDATE METADATA
#-------------------------------------------------------------------------------
echo ""
echo "Validating sample metadata..."

qiime metadata tabulate \
    --m-input-file "${METADATA_DIR}/sample_metadata.tsv" \
    --o-visualization "${QIIME2_DIR}/imported/metadata-summary.qzv"

check_status "Metadata validation"

# conda deactivate

#-------------------------------------------------------------------------------
# SUMMARY
#-------------------------------------------------------------------------------
echo ""
echo "=============================================="
echo "QIIME2 Import Complete!"
echo "=============================================="
echo "End time: $(date)"
echo ""
echo "Output files:"
echo "  Manifest: ${MANIFEST_FILE}"
echo "  Imported sequences: ${OUTPUT_QZA}"
echo "  Demux summary: ${DEMUX_VIZ}"
echo "  Exported data: ${DEMUX_EXPORT}/"
echo ""
echo "IMPORTANT: Open the demux summary to determine DADA2 parameters:"
echo "  View with: qiime tools view ${DEMUX_VIZ}"
echo "  Or open: ${DEMUX_EXPORT}/per-sample-fastq-counts.tsv"
echo ""
echo "Check the quality plots to determine:"
echo "  - Forward truncation length (--p-trunc-len-f)"
echo "  - Reverse truncation length (--p-trunc-len-r)"
echo ""
echo "Update config.sh with these values before running DADA2:"
echo "  TRUNC_LEN_F=<value>"
echo "  TRUNC_LEN_R=<value>"
echo ""
echo "Next step: bash ${SCRIPT_DIR}/05_dada2_denoise.sh"

duration=$SECONDS
echo "$((duration / 60)) minutes $((duration % 60)) seconds elapsed."