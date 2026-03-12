#!/bin/bash
#===============================================================================
# 00_SETUP_ENVIRONMENT.SH - Environment and Database Setup
#===============================================================================
# This script sets up the computational environment for the 16S analysis pipeline
# Run this script FIRST before any other scripts
#===============================================================================

set -e  # Exit on error

# Source configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

echo "=============================================="
echo "16S Microbiome Pipeline - Environment Setup"
echo "=============================================="
echo "Start time: $(date)"
echo ""

#-------------------------------------------------------------------------------
# CREATE DIRECTORY STRUCTURE
#-------------------------------------------------------------------------------
echo "[1/5] Creating directory structure..."
create_directories

#-------------------------------------------------------------------------------
# CHECK/INSTALL CONDA ENVIRONMENTS
#-------------------------------------------------------------------------------
echo ""
echo "[2/5] Setting up conda environments..."

# Check if conda is available
if ! command -v conda &> /dev/null; then
    echo "ERROR: conda not found. Please install Miniconda or Anaconda first."
    echo "Visit: https://docs.conda.io/en/latest/miniconda.html"
    exit 1
fi

# Initialize conda for bash if needed
eval "$(conda shell.bash hook)"

# Function to create or update conda env
setup_conda_env() {
    local env_name=$1
    local env_file=$2
    
    if conda env list | grep -q "^${env_name} "; then
        echo "  Environment '${env_name}' already exists. Skipping..."
    else
        echo "  Creating environment '${env_name}' from ${env_file}..."
        if [ -f "${env_file}" ]; then
            conda env create -f "${env_file}" || {
                echo "  Warning: Could not create from yml. Will use alternative method."
            }
        fi
    fi
}

# QIIME2 environment - recommend official installer
echo ""
echo "  Checking QIIME2 environment..."
if conda env list | grep -q "^${QIIME2_ENV} "; then
    echo "  QIIME2 environment '${QIIME2_ENV}' exists."
else
    echo "  QIIME2 environment not found."
    echo "  For best results, install QIIME2 using the official installer:"
    echo "  wget https://raw.githubusercontent.com/qiime2/distributions/refs/heads/dev/2026.1/amplicon/released/qiime2-amplicon-ubuntu-latest-conda.yml"
    echo "  conda env create -n ${QIIME2_ENV} --file qiime2-amplicon-ubuntu-latest-conda.yml"
    echo ""
    echo "  Alternatively, attempting to create directly from URL..."
    
    # Try creating directly from URL or fall back to local file if setup_conda_env was structured differently
    # But here let's just run the command directly since we have the specific user request
    conda env create \
      --name "${QIIME2_ENV}" \
      --file "https://raw.githubusercontent.com/qiime2/distributions/refs/heads/dev/2026.1/amplicon/released/qiime2-amplicon-ubuntu-latest-conda.yml" || {
        echo "  Warning: Could not create from URL. Please install manually."
    }
fi

# Install FastQC and MultiQC in QIIME2 env if not present
echo "  Ensuring FastQC and MultiQC are installed..."
conda activate "${QIIME2_ENV}" 2>/dev/null || true
if command -v fastqc &> /dev/null; then
    echo "    FastQC: $(fastqc --version 2>&1 | head -1)"
else
    echo "    Installing FastQC..."
    conda install -c bioconda fastqc -y 2>/dev/null || echo "    Note: Install FastQC manually if needed"
fi
if command -v multiqc &> /dev/null; then
    echo "    MultiQC: $(multiqc --version 2>&1)"
else
    echo "    Installing MultiQC..."
    conda install -c bioconda multiqc -y 2>/dev/null || echo "    Note: Install MultiQC manually if needed"
fi
conda deactivate 2>/dev/null || true

# R environment
echo ""
echo "  Checking R environment..."
setup_conda_env "${R_ENV}" "${PROJECT_DIR}/envs/r-microbiome-env.yml"

# PICRUSt2 environment  
echo ""
echo "  Checking PICRUSt2 environment..."
setup_conda_env "${PICRUST2_ENV}" "${PROJECT_DIR}/envs/picrust2-env.yml"

#-------------------------------------------------------------------------------
# INSTALL ADDITIONAL R PACKAGES
#-------------------------------------------------------------------------------
echo ""
echo "[3/5] Installing additional R packages..."

# Create R package installation script
cat > "${SCRIPT_DIR}/install_r_packages.R" << 'RSCRIPT'
# Install required R packages not available via conda

# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Install BiocManager if needed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

# Function to install if not present
install_if_missing <- function(pkg, bioc = FALSE) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        message(paste("Installing", pkg, "..."))
        if (bioc) {
            BiocManager::install(pkg, ask = FALSE, update = FALSE)
        } else {
            install.packages(pkg)
        }
    } else {
        message(paste(pkg, "already installed"))
    }
}

# Bioconductor packages
install_if_missing("ANCOMBC", bioc = TRUE)
install_if_missing("ALDEx2", bioc = TRUE)
install_if_missing("phyloseq", bioc = TRUE)
install_if_missing("microbiome", bioc = TRUE)
install_if_missing("ComplexHeatmap", bioc = TRUE)
install_if_missing("DESeq2", bioc = TRUE)

# GitHub packages
if (!requireNamespace("qiime2R", quietly = TRUE)) {
    message("Installing qiime2R from GitHub...")
    if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
    devtools::install_github("jbisanz/qiime2R")
}

if (!requireNamespace("microViz", quietly = TRUE)) {
    message("Installing microViz from GitHub...")
    if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
    remotes::install_github("david-barnett/microViz")
}

message("\nR package installation complete!")
RSCRIPT

echo "  R package installation script created at: ${SCRIPT_DIR}/install_r_packages.R"
echo "  Run this after activating R environment:"
echo "    conda activate ${R_ENV}"
echo "    Rscript ${SCRIPT_DIR}/install_r_packages.R"

#-------------------------------------------------------------------------------
# DOWNLOAD TAXONOMY DATABASES
#-------------------------------------------------------------------------------
echo ""
echo "[4/5] Setting up taxonomy databases..."

# SILVA database
echo ""
echo "  SILVA ${SILVA_VERSION} database..."
if [ -f "${SILVA_CLASSIFIER}" ]; then
    echo "    SILVA classifier already exists."
else
    echo "    SILVA classifier will be automatically downloaded and trained in step 06."
    echo "    Source URLs configured in config.sh"
fi

# GTDB database
echo ""
echo "  GTDB ${GTDB_VERSION} database..."
if [ -f "${GTDB_CLASSIFIER}" ]; then
    echo "    GTDB classifier already exists."
else
    echo "    Note: GTDB classifier needs to be trained or downloaded separately."
    echo "    Options:"
    echo "    1. Use RESCRIPt to train classifier (recommended):"
    echo "       qiime rescript get-gtdb-data --p-version '214.1' ..."
    echo "    2. Download pre-trained classifier if available"
    echo ""
    echo "    For now, creating placeholder. Run training script when ready."
    touch "${DATABASE_DIR}/gtdb/DOWNLOAD_GTDB_CLASSIFIER.txt"
    echo "Download GTDB classifier or train using RESCRIPt" > "${DATABASE_DIR}/gtdb/DOWNLOAD_GTDB_CLASSIFIER.txt"
fi

#-------------------------------------------------------------------------------
# VERIFY SETUP
#-------------------------------------------------------------------------------
echo ""
echo "[5/5] Verifying setup..."

echo ""
echo "Directory structure:"
echo "  Results: $(ls -d ${RESULTS_DIR} 2>/dev/null && echo 'OK' || echo 'MISSING')"
echo "  QC Reports: $(ls -d ${QC_DIR} 2>/dev/null && echo 'OK' || echo 'MISSING')"
echo "  QIIME2: $(ls -d ${QIIME2_DIR} 2>/dev/null && echo 'OK' || echo 'MISSING')"
echo "  Figures: $(ls -d ${FIGURES_DIR} 2>/dev/null && echo 'OK' || echo 'MISSING')"

echo ""
echo "Conda environments:"
conda env list | grep -E "(${QIIME2_ENV}|${R_ENV}|${PICRUST2_ENV})" || echo "  Some environments may need manual setup"

echo ""
echo "Databases:"
echo "  SILVA: $([ -f "${SILVA_CLASSIFIER}" ] && echo 'OK' || echo 'Will be trained automatically')"
echo "  GTDB: $([ -f "${GTDB_CLASSIFIER}" ] && echo 'OK' || echo 'Not downloaded')"

echo ""
echo "=============================================="
echo "Setup complete!"
echo "End time: $(date)"
echo "=============================================="
echo ""
echo "Next steps:"
echo "1. If any conda environments failed, install manually"
echo "2. Activate R environment and run: Rscript ${SCRIPT_DIR}/install_r_packages.R"
echo "3. Download/train GTDB classifier if needed"
echo "4. Run: bash ${SCRIPT_DIR}/01_fastqc_pretrim.sh"
