#!/bin/bash
#===============================================================================
# 10_PICRUST2.SH - Functional Prediction with PICRUSt2
#===============================================================================
# Predict metagenomic functional content from 16S data
#===============================================================================
SECONDS=0
set -e

# Source configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

echo "=============================================="
echo "PICRUSt2 Functional Prediction"
echo "=============================================="
echo "Start time: $(date)"
echo ""

#Check if PICRUSt2 environment is active
if [[ "${CONDA_DEFAULT_ENV}" != "${PICRUST2_ENV}" ]]; then
    echo "ERROR: PICRUSt2 environment is not active."
    echo "Expected environment: ${PICRUST2_ENV}"
    echo "Current environment: ${CONDA_DEFAULT_ENV:-None}"
    echo ""
    echo "Please activate the environment before running this script:"
    echo "  conda activate ${PICRUST2_ENV}"
    exit 1
fi

# Verify PICRUSt2 is installed
picrust2_pipeline.py --version || {
    echo "ERROR: PICRUSt2 not found. Please install it first."
    echo "conda create -n picrust2 -c bioconda -c conda-forge picrust2"
    exit 1
}

OUTPUT_DIR="${PICRUST2_DIR}"
mkdir -p "${OUTPUT_DIR}"

#-------------------------------------------------------------------------------
# PREPARE INPUT FILES
#-------------------------------------------------------------------------------
echo "Preparing input files..."

# Input files from DADA2
ASV_SEQS="${TABLES_DIR}/asv_sequences.fasta"
ASV_TABLE="${TABLES_DIR}/asv_counts.tsv"

# Check inputs exist
if [ ! -f "${ASV_SEQS}" ]; then
    echo "ERROR: ASV sequences not found: ${ASV_SEQS}"
    echo "Run DADA2 pipeline first."
    exit 1
fi

if [ ! -f "${ASV_TABLE}" ]; then
    echo "ERROR: ASV table not found: ${ASV_TABLE}"
    echo "Run DADA2 pipeline first."
    exit 1
fi

# Safely format the header for BIOM conversion without losing sample IDs.
# The original `grep -v "^#"` stripped ALL lines starting with #, including
# the header row that contains sample names (e.g. "#OTU ID  C1  C2  C3").
# This awk command preserves sample names while ensuring the first column
# is named "#OTU ID" as required by biom convert.
awk 'NR==1 {
    if ($1 !~ /^#/) { $1="#OTU ID" }
    print
}
NR>1 {print}' OFS='\t' "${ASV_TABLE}" > "${OUTPUT_DIR}/asv_counts_clean.tsv"

# Convert to BIOM
biom convert \
    -i "${OUTPUT_DIR}/asv_counts_clean.tsv" \
    -o "${OUTPUT_DIR}/asv_table.biom" \
    --table-type="OTU table" \
    --to-hdf5

check_status "BIOM conversion"

#-------------------------------------------------------------------------------
# RUN PICRUST2 FULL PIPELINE
#-------------------------------------------------------------------------------
echo ""
echo "Running PICRUSt2 pipeline..."
echo "This typically takes 30 minutes to 2 hours depending on data size."
echo ""

picrust2_pipeline.py \
    -s "${ASV_SEQS}" \
    -i "${OUTPUT_DIR}/asv_table.biom" \
    -o "${OUTPUT_DIR}/picrust2_out" \
    -p ${THREADS} \
    --in_traits EC,KO \
    --stratified \
    --coverage \
    --verbose

check_status "PICRUSt2 pipeline"

#-------------------------------------------------------------------------------
# ADD PATHWAY DESCRIPTIONS
#-------------------------------------------------------------------------------
echo ""
echo "Adding descriptions to pathway predictions..."

# Add EC descriptions
if [ -f "${OUTPUT_DIR}/picrust2_out/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz" ]; then
    gunzip -k "${OUTPUT_DIR}/picrust2_out/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz" 2>/dev/null || true
    
    add_descriptions.py \
        -i "${OUTPUT_DIR}/picrust2_out/EC_metagenome_out/pred_metagenome_unstrat.tsv" \
        -m EC \
        -o "${OUTPUT_DIR}/EC_metagenome_with_descriptions.tsv"
fi

# Add KO descriptions
if [ -f "${OUTPUT_DIR}/picrust2_out/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz" ]; then
    gunzip -k "${OUTPUT_DIR}/picrust2_out/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz" 2>/dev/null || true
    
    add_descriptions.py \
        -i "${OUTPUT_DIR}/picrust2_out/KO_metagenome_out/pred_metagenome_unstrat.tsv" \
        -m KO \
        -o "${OUTPUT_DIR}/KO_metagenome_with_descriptions.tsv"
fi

# Add pathway descriptions
if [ -f "${OUTPUT_DIR}/picrust2_out/pathways_out/path_abun_unstrat.tsv.gz" ]; then
    gunzip -k "${OUTPUT_DIR}/picrust2_out/pathways_out/path_abun_unstrat.tsv.gz" 2>/dev/null || true
    
    add_descriptions.py \
        -i "${OUTPUT_DIR}/picrust2_out/pathways_out/path_abun_unstrat.tsv" \
        -m METACYC \
        -o "${OUTPUT_DIR}/pathways_with_descriptions.tsv"
fi

check_status "Description addition"

#-------------------------------------------------------------------------------
# NSTI QUALITY CHECK
#-------------------------------------------------------------------------------
echo ""
echo "Checking NSTI (Nearest Sequenced Taxon Index) values..."

NSTI_FILE="${OUTPUT_DIR}/picrust2_out/out.tre"
if [ -f "${OUTPUT_DIR}/picrust2_out/EC_metagenome_out/weighted_nsti.tsv.gz" ]; then
    gunzip -k "${OUTPUT_DIR}/picrust2_out/EC_metagenome_out/weighted_nsti.tsv.gz" 2>/dev/null || true
    NSTI_FILE="${OUTPUT_DIR}/picrust2_out/EC_metagenome_out/weighted_nsti.tsv"
fi

if [ -f "${NSTI_FILE}" ]; then
    echo "NSTI Statistics:"
    # Calculate NSTI summary
    awk -F'\t' 'NR>1 {sum+=$2; count++; if($2>max)max=$2; if(min==""||$2<min)min=$2} 
         END {print "  Mean NSTI:", sum/count; print "  Min NSTI:", min; print "  Max NSTI:", max}' \
        "${NSTI_FILE}"
    
    echo ""
    echo "Note: NSTI > 2 indicates poor prediction accuracy for those samples."
    echo "Samples with high NSTI should be interpreted with caution."
fi

#-------------------------------------------------------------------------------
# ORGANIZE OUTPUTS
#-------------------------------------------------------------------------------
echo ""
echo "Organizing output files..."

# Copy key files to main results directory
cp "${OUTPUT_DIR}/EC_metagenome_with_descriptions.tsv" "${TABLES_DIR}/" 2>/dev/null || true
cp "${OUTPUT_DIR}/KO_metagenome_with_descriptions.tsv" "${TABLES_DIR}/" 2>/dev/null || true
cp "${OUTPUT_DIR}/pathways_with_descriptions.tsv" "${TABLES_DIR}/" 2>/dev/null || true

# Copy stratified outputs if generated
if [ -d "${OUTPUT_DIR}/picrust2_out/pathways_out" ]; then
    mkdir -p "${TABLES_DIR}/picrust2_stratified"
    cp ${OUTPUT_DIR}/picrust2_out/pathways_out/*strat*.tsv* \
       "${TABLES_DIR}/picrust2_stratified/" 2>/dev/null || true
fi

check_status "Output organization"

#-------------------------------------------------------------------------------
# SUMMARY
#-------------------------------------------------------------------------------
echo ""
echo "=============================================="
echo "PICRUSt2 Analysis Complete!"
echo "=============================================="
echo "End time: $(date)"
echo ""
echo "Output files:"
echo "  EC predictions: ${OUTPUT_DIR}/EC_metagenome_with_descriptions.tsv"
echo "  KO predictions: ${OUTPUT_DIR}/KO_metagenome_with_descriptions.tsv"
echo "  Pathway predictions: ${OUTPUT_DIR}/pathways_with_descriptions.tsv"
echo "  Full output: ${OUTPUT_DIR}/picrust2_out/"
echo ""
echo "Key outputs for downstream analysis:"
echo "  - EC_metagenome: Enzyme Commission numbers"
echo "  - KO_metagenome: KEGG Orthology groups"
echo "  - pathways: MetaCyc pathway abundances"
echo ""
echo "Next step: Run R analysis scripts in ${R_ANALYSIS_DIR}/"

duration=$SECONDS
echo "Duration: $(($duration / 60)) minutes and $(($duration % 60)) seconds"