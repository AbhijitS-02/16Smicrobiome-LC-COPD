#!/bin/bash
#===============================================================================
# CONFIG.SH - Central Configuration for 16S Microbiome Analysis Pipeline
#===============================================================================
# Project: LC-COPD Microbiome Study
# Author: Abhijit
# Date: 2026-02-08
#===============================================================================

#-------------------------------------------------------------------------------
# DIRECTORY PATHS
#-------------------------------------------------------------------------------
export PROJECT_DIR="/home/abhijit/LC_COPD_microbiome"
export RAW_DATA_DIR="${PROJECT_DIR}/raw_data"
export SCRIPTS_DIR="${PROJECT_DIR}/scripts"
export RESULTS_DIR="${PROJECT_DIR}/results"
export FIGURES_DIR="${PROJECT_DIR}/figures"
export METADATA_DIR="${PROJECT_DIR}/metadata"
export DATABASE_DIR="${PROJECT_DIR}/databases"
export LOGS_DIR="${PROJECT_DIR}/logs"
export R_ANALYSIS_DIR="${PROJECT_DIR}/R_analysis"

# Results subdirectories
export QC_DIR="${RESULTS_DIR}/qc_reports"
export TRIMMED_DIR="${RESULTS_DIR}/trimmed_reads"
export QIIME2_DIR="${RESULTS_DIR}/qiime2"
export PICRUST2_DIR="${RESULTS_DIR}/picrust2"
export TABLES_DIR="${RESULTS_DIR}/tables"

#-------------------------------------------------------------------------------
# COMPUTATIONAL RESOURCES
#-------------------------------------------------------------------------------
export THREADS=6          # Ryzen 5 8645HS (6 cores, using 6 threads for stability)
export MAX_MEMORY="10G"   # ~11GB available, leaving buffer

#-------------------------------------------------------------------------------
# PRIMER SEQUENCES (V3-V4 Region: 341F with 805R/806R)
#-------------------------------------------------------------------------------
# Forward primer: 341F (5'-CCTACGGGNGGCWGCAG-3') - 17bp
# Reverse primers: Data shows mixed usage:
#   - 805R (5'-GACTACHVGGGTATCTAATCC-3') - 21bp (39.3% detection)
#   - 806R (5'-GGACTACHVGGGTWTCTAAT-3') - 20bp (89.0% detection)
# Both will be trimmed in the cutadapt step
export FORWARD_PRIMER="CCTACGGGNGGCWGCAG"
export REVERSE_PRIMER_805R="GACTACHVGGGTATCTAATCC"
export REVERSE_PRIMER_806R="GGACTACHVGGGTWTCTAAT"
# Primary reverse primer (for scripts that only use one)
export REVERSE_PRIMER="${REVERSE_PRIMER_806R}"

# Trimming parameters based on primer detection analysis
# 341F detected in R1 (96.1%), reverse primers partially detected
export TRIM_LEFT_F=0     # Trim 341F primer from forward reads (can be 0 since primers already trimmed using cutadapt)
export TRIM_LEFT_R=0      # Reverse primers appear pre-trimmed

#-------------------------------------------------------------------------------
# QUALITY FILTERING PARAMETERS
#-------------------------------------------------------------------------------
export MIN_QUALITY=20     # Minimum Phred quality score
export MIN_LENGTH=200     # Minimum read length after trimming
export MAX_N=0            # Maximum N bases allowed

#-------------------------------------------------------------------------------
# DADA2 PARAMETERS
#-------------------------------------------------------------------------------
# These will be adjusted after reviewing FastQC quality profiles
# TRUNC_LEN and MAX_EE parameters are set after running optimize_dada2_params.R
# A grid-search method was used with variable truncation length and maximum expected errors
# trunc_f_opts <- 260, 270, 275; trunc_r_opts <- 210, 220; max_ee_opts <- (2, 4), (3, 4)
# Error learning curves were plotted and selection of best parameters was based on :
# 1. Convergence of error rate to a stable value
# 2. Overlap of >=20bp for an amplicon length of ~460bp (V3-V4 region)
# 3. Modest max expected errors (slightly relaxed for reverse reads)

export TRUNC_LEN_F=275      # To be set after QC (0 = no truncation)
export TRUNC_LEN_R=220      # To be set after QC (0 = no truncation)
export MAX_EE_F=2         # Maximum expected errors (forward)
export MAX_EE_R=4         # Maximum expected errors (reverse)
export TRUNC_Q=2          # Truncate at first base with quality <= this
export CHIMERA_METHOD="consensus"
export MIN_FOLD_PARENT=1  # For chimera detection

#-------------------------------------------------------------------------------
# TAXONOMY DATABASES
#-------------------------------------------------------------------------------
export SILVA_VERSION="138-99-515-806"
export SILVA_CLASSIFIER="${DATABASE_DIR}/silva/silva-${SILVA_VERSION}-nb-classifier.qza"
export SILVA_SEQS_URL="https://data.qiime2.org/2024.5/common/silva-138-99-seqs.qza"
export SILVA_TAX_URL="https://data.qiime2.org/2024.5/common/silva-138-99-tax.qza"

export GTDB_VERSION="r214"
export GTDB_CLASSIFIER="${DATABASE_DIR}/gtdb/gtdb-${GTDB_VERSION}-ssu-nb-classifier.qza"

#-------------------------------------------------------------------------------
# DIVERSITY ANALYSIS PARAMETERS
#-------------------------------------------------------------------------------
export RAREFACTION_DEPTH=0  # To be set after feature table summary (0 = auto)
export ALPHA_METRICS="observed_features shannon simpson faith_pd pielou_e chao1"
export BETA_METRICS="bray_curtis jaccard weighted_unifrac unweighted_unifrac"

#-------------------------------------------------------------------------------
# RANDOM SEED FOR REPRODUCIBILITY
#-------------------------------------------------------------------------------
export RANDOM_SEED=42

#-------------------------------------------------------------------------------
# SAMPLE GROUPS
#-------------------------------------------------------------------------------
export GROUP_CONTROL="Control"
export GROUP_LC_COPD="LC_COPD"
export GROUP_LC_WOCOPD="LC_woCOPD"

#-------------------------------------------------------------------------------
# CONDA ENVIRONMENTS
#-------------------------------------------------------------------------------
export QIIME2_ENV="qiime2-amplicon-2026.1"
export R_ENV="microbiome-r"
export PICRUST2_ENV="picrust2"

#-------------------------------------------------------------------------------
# FILE PATTERNS
#-------------------------------------------------------------------------------
export R1_PATTERN="_R1.fastq.gz"
export R2_PATTERN="_R2.fastq.gz"

#-------------------------------------------------------------------------------
# LOGGING CONFIGURATION
#-------------------------------------------------------------------------------
export LOG_DATE=$(date +"%Y%m%d_%H%M%S")

# Function to log messages
log_message() {
    local level=$1
    local message=$2
    local timestamp=$(date +"%Y-%m-%d %H:%M:%S")
    echo "[${timestamp}] [${level}] ${message}" | tee -a "${LOGS_DIR}/pipeline_${LOG_DATE}.log"
}

# Function to check if command succeeded
check_status() {
    if [ $? -eq 0 ]; then
        log_message "INFO" "$1 completed successfully"
    else
        log_message "ERROR" "$1 failed"
        exit 1
    fi
}

#-------------------------------------------------------------------------------
# CREATE DIRECTORY STRUCTURE
#-------------------------------------------------------------------------------
create_directories() {
    # Create logs directory first so log_message can write to log file
    mkdir -p "${LOGS_DIR}"
    
    log_message "INFO" "Creating directory structure..."
    mkdir -p "${RESULTS_DIR}"
    mkdir -p "${QC_DIR}/pretrim_fastqc"
    mkdir -p "${QC_DIR}/posttrim_fastqc"
    mkdir -p "${QC_DIR}/multiqc"
    mkdir -p "${TRIMMED_DIR}"
    mkdir -p "${QIIME2_DIR}/imported"
    mkdir -p "${QIIME2_DIR}/dada2_output"
    mkdir -p "${QIIME2_DIR}/taxonomy_silva"
    mkdir -p "${QIIME2_DIR}/taxonomy_gtdb"
    mkdir -p "${QIIME2_DIR}/phylogeny"
    mkdir -p "${QIIME2_DIR}/diversity"
    mkdir -p "${PICRUST2_DIR}"
    mkdir -p "${TABLES_DIR}"
    mkdir -p "${FIGURES_DIR}/qc"
    mkdir -p "${FIGURES_DIR}/alpha_diversity"
    mkdir -p "${FIGURES_DIR}/beta_diversity"
    mkdir -p "${FIGURES_DIR}/composition"
    mkdir -p "${FIGURES_DIR}/differential_abundance"
    mkdir -p "${FIGURES_DIR}/functional"
    mkdir -p "${METADATA_DIR}"
    mkdir -p "${DATABASE_DIR}/silva"
    mkdir -p "${DATABASE_DIR}/gtdb"
    mkdir -p "${R_ANALYSIS_DIR}/functions"
    
    check_status "Directory structure creation"
}

echo "Configuration loaded from: ${BASH_SOURCE[0]}"
echo "Project directory: ${PROJECT_DIR}"
echo "Using ${THREADS} threads with ${MAX_MEMORY} memory limit"
