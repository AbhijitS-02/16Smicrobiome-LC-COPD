#!/bin/bash
#===============================================================================
# 02_TRIM_PRIMERS.SH - Primer Trimming and Quality Filtering
#===============================================================================
# Remove primers using cutadapt and perform quality filtering
# V3-V4 region: 341F forward with 805R/806R reverse primers
# Note: Data shows mixed reverse primers (805R: 39.3%, 806R: 89.0% detection)
#===============================================================================
SECONDS=0
set -e

# Source configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

echo "=============================================="
echo "Primer Trimming and Quality Filtering"
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
# CONFIGURATION
#-------------------------------------------------------------------------------
echo "Trimming Configuration:"
echo "  Forward primer (341F): ${FORWARD_PRIMER} (${#FORWARD_PRIMER} bp)"
echo "  Reverse primer (805R): ${REVERSE_PRIMER_805R} (${#REVERSE_PRIMER_805R} bp)"
echo "  Reverse primer (806R): ${REVERSE_PRIMER_806R} (${#REVERSE_PRIMER_806R} bp)"
echo "  Trim left F: ${TRIM_LEFT_F}"
echo "  Trim left R: ${TRIM_LEFT_R}"
echo "  Min quality: ${MIN_QUALITY}"
echo "  Min length: ${MIN_LENGTH}"
echo ""

# Create output directory
mkdir -p "${TRIMMED_DIR}"
mkdir -p "${LOGS_DIR}/trimming"

#-------------------------------------------------------------------------------
# NEXTERA TRANSPOSASE ADAPTER SEQUENCES
#-------------------------------------------------------------------------------
# These adapters were detected in MultiQC adapter content analysis
# They appear at 3' end of reads and need to be removed

# Nextera Read 1 transposase: TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG
# Nextera Read 2 transposase: GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG
# Also the reverse complements for reads that sequence through adapters

NEXTERA_R1="TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"
NEXTERA_R2="GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG"
# Reverse complements
NEXTERA_R1_RC="CTGTCTCTTATACACATCTGACGCTGCCGACGA"
NEXTERA_R2_RC="CTGTCTCTTATACACATCTCCGAGCCCACGAGAC"

echo "Adapter sequences to remove:"
echo "  Forward primer (341F): ${FORWARD_PRIMER}"
echo "  Reverse primer (805R): ${REVERSE_PRIMER_805R}"
echo "  Reverse primer (806R): ${REVERSE_PRIMER_806R}"
echo "  Nextera R1: ${NEXTERA_R1}"
echo "  Nextera R2: ${NEXTERA_R2}"
echo ""

# Compute reverse complements for 3' adapter trimming
# 805R RC: GGATTAGATACCCTBHGTAGTC (approximate, using the standard RC)
# 806R RC: ATTAGAWACCCBHGTAGTCC (approximate)
REVERSE_PRIMER_805R_RC="GGATTAGATACCCBDGTAGTC"
REVERSE_PRIMER_806R_RC="ATTAGAWACCCBDGTAGTCC"

#-------------------------------------------------------------------------------
# TRIM PRIMERS AND ADAPTERS WITH CUTADAPT
#-------------------------------------------------------------------------------
echo "Running cutadapt for primer and adapter removal..."
echo ""

# Get list of samples
SAMPLES=$(ls ${RAW_DATA_DIR}/*_R1.fastq.gz | xargs -n1 basename | sed 's/_R1.fastq.gz//')

# Create manifest for tracking
MANIFEST="${TRIMMED_DIR}/trimming_manifest.txt"
echo -e "sample-id\tinput_R1\tinput_R2\toutput_R1\toutput_R2" > "${MANIFEST}"

# Process each sample
TOTAL=$(echo "${SAMPLES}" | wc -w)
COUNT=0

for SAMPLE in ${SAMPLES}; do
    COUNT=$((COUNT + 1))
    echo "[${COUNT}/${TOTAL}] Processing ${SAMPLE}..."
    
    INPUT_R1="${RAW_DATA_DIR}/${SAMPLE}_R1.fastq.gz"
    INPUT_R2="${RAW_DATA_DIR}/${SAMPLE}_R2.fastq.gz"
    OUTPUT_R1="${TRIMMED_DIR}/${SAMPLE}_R1_trimmed.fastq.gz"
    OUTPUT_R2="${TRIMMED_DIR}/${SAMPLE}_R2_trimmed.fastq.gz"
    LOG_FILE="${LOGS_DIR}/trimming/${SAMPLE}_cutadapt.log"
    
    # Run cutadapt
    # -g: 5' adapter on R1 (forward primer - non-anchored)
    # -G: 5' adapter on R2 (reverse primer - non-anchored)
    # -a: 3' adapters on R1 (reverse primer RC + Nextera adapters)
    # -A: 3' adapters on R2 (forward primer RC + Nextera adapters)
    # --minimum-length: discard reads shorter than this
    # -q: quality trimming threshold (from both ends)
    # --trim-n: remove N bases from ends
    # -n 2: remove up to 2 adapters from each read (primer + Nextera)
    
    # Trim both 805R and 806R reverse primers (data shows mixed usage)
    # -g: 5' adapter on R1 (forward primer)
    # -G: 5' adapters on R2 (both reverse primers - cutadapt tries each)
    # -a: 3' adapters on R1 (reverse primer RCs)
    # -A: 3' adapters on R2 (forward primer RC)
    cutadapt \
        -g "${FORWARD_PRIMER}" \
        -G "${REVERSE_PRIMER_805R}" \
        -G "${REVERSE_PRIMER_806R}" \
        -a "${REVERSE_PRIMER_805R_RC}" \
        -a "${REVERSE_PRIMER_806R_RC}" \
        -a "${NEXTERA_R1_RC}" \
        -a "${NEXTERA_R2_RC}" \
        -A "${FORWARD_PRIMER_RC:-CTGCWGCCNCCCGTAGG}" \
        -A "${NEXTERA_R1_RC}" \
        -A "${NEXTERA_R2_RC}" \
        -n 3 \
        --minimum-length ${MIN_LENGTH} \
        -q ${MIN_QUALITY},${MIN_QUALITY} \
        --trim-n \
        --cores ${THREADS} \
        -o "${OUTPUT_R1}" \
        -p "${OUTPUT_R2}" \
        "${INPUT_R1}" \
        "${INPUT_R2}" \
        > "${LOG_FILE}" 2>&1
    
    # Add to manifest
    echo -e "${SAMPLE}\t${INPUT_R1}\t${INPUT_R2}\t${OUTPUT_R1}\t${OUTPUT_R2}" >> "${MANIFEST}"
    
done

check_status "Cutadapt trimming"

#-------------------------------------------------------------------------------
# GENERATE TRIMMING STATISTICS
#-------------------------------------------------------------------------------
echo ""
echo "Generating trimming statistics..."

STATS_FILE="${TRIMMED_DIR}/trimming_stats.txt"

echo "Trimming Statistics Summary" > "${STATS_FILE}"
echo "Generated: $(date)" >> "${STATS_FILE}"
echo "==============================" >> "${STATS_FILE}"
echo "" >> "${STATS_FILE}"

# Parse cutadapt logs for statistics
echo -e "Sample\tTotal_Reads\tPassing_Reads\tPercent_Passing" >> "${STATS_FILE}"

for LOG in ${LOGS_DIR}/trimming/*_cutadapt.log; do
    SAMPLE=$(basename "${LOG}" _cutadapt.log)
    
    # Extract statistics from cutadapt log
    TOTAL_READS=$(grep "Total read pairs processed:" "${LOG}" | awk '{gsub(/,/,""); print $NF}')
    PASSING_READS=$(grep "Pairs written (passing filters):" "${LOG}" | awk '{gsub(/,/,""); print $5}')
    
    if [ -n "${TOTAL_READS}" ] && [ -n "${PASSING_READS}" ]; then
        PERCENT=$(echo "scale=2; ${PASSING_READS}/${TOTAL_READS}*100" | bc)
        echo -e "${SAMPLE}\t${TOTAL_READS}\t${PASSING_READS}\t${PERCENT}%" >> "${STATS_FILE}"
    fi
done

echo "" >> "${STATS_FILE}"

# Calculate overall statistics
echo "Overall Summary:" >> "${STATS_FILE}"
TOTAL_IN=$(cat ${LOGS_DIR}/trimming/*.log | grep "Total read pairs processed:" | awk '{gsub(/,/,""); sum+=$NF} END{print sum}')
TOTAL_OUT=$(cat ${LOGS_DIR}/trimming/*.log | grep "Pairs written (passing filters):" | awk '{gsub(/,/,""); sum+=$5} END{print sum}')
OVERALL_PERCENT=$(echo "scale=2; ${TOTAL_OUT}/${TOTAL_IN}*100" | bc)

echo "  Total input read pairs: ${TOTAL_IN}" >> "${STATS_FILE}"
echo "  Total output read pairs: ${TOTAL_OUT}" >> "${STATS_FILE}"
echo "  Overall pass rate: ${OVERALL_PERCENT}%" >> "${STATS_FILE}"

cat "${STATS_FILE}"

# conda deactivate

#-------------------------------------------------------------------------------
# SUMMARY
#-------------------------------------------------------------------------------
echo ""
echo "=============================================="
echo "Primer Trimming Complete!"
echo "=============================================="
echo "End time: $(date)"
echo ""
echo "Output files:"
echo "  Trimmed reads: ${TRIMMED_DIR}/"
echo "  Manifest: ${MANIFEST}"
echo "  Statistics: ${STATS_FILE}"
echo "  Logs: ${LOGS_DIR}/trimming/"
echo ""
echo "Next step: bash ${SCRIPT_DIR}/03_fastqc_posttrim.sh"

duration=$SECONDS
echo "$((duration / 60)) minutes $((duration % 60)) seconds elapsed."