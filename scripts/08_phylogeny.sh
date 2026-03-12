#!/bin/bash
#===============================================================================
# 08_PHYLOGENY.SH - Phylogenetic Tree Construction
#===============================================================================
# Build a phylogenetic tree from ASV sequences for diversity metrics
#===============================================================================
SECONDS=0
set -e

# Source configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

echo "=============================================="
echo "Phylogenetic Tree Construction"
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

INPUT_SEQS="${QIIME2_DIR}/dada2_output/rep-seqs.qza"
OUTPUT_DIR="${QIIME2_DIR}/phylogeny"

mkdir -p "${OUTPUT_DIR}"

#-------------------------------------------------------------------------------
# OPTION 1: DE NOVO PHYLOGENY (DEFAULT)
#-------------------------------------------------------------------------------
echo "Building de novo phylogenetic tree..."
echo ""

# Multiple sequence alignment with MAFFT
echo "Step 1/4: Aligning sequences with MAFFT..."
qiime alignment mafft \
    --i-sequences "${INPUT_SEQS}" \
    --p-n-threads ${THREADS} \
    --o-alignment "${OUTPUT_DIR}/aligned-rep-seqs.qza"

check_status "MAFFT alignment"

# Mask highly variable regions
echo "Step 2/4: Masking highly variable positions..."
qiime alignment mask \
    --i-alignment "${OUTPUT_DIR}/aligned-rep-seqs.qza" \
    --o-masked-alignment "${OUTPUT_DIR}/masked-aligned-rep-seqs.qza"

check_status "Alignment masking"

# Build tree with FastTree
echo "Step 3/4: Building tree with FastTree..."
qiime phylogeny fasttree \
    --i-alignment "${OUTPUT_DIR}/masked-aligned-rep-seqs.qza" \
    --p-n-threads ${THREADS} \
    --o-tree "${OUTPUT_DIR}/unrooted-tree.qza"

check_status "FastTree construction"

# Root tree at midpoint
echo "Step 4/4: Rooting tree at midpoint..."
qiime phylogeny midpoint-root \
    --i-tree "${OUTPUT_DIR}/unrooted-tree.qza" \
    --o-rooted-tree "${OUTPUT_DIR}/rooted-tree.qza"

check_status "Tree rooting"

#-------------------------------------------------------------------------------
# ALTERNATIVE: FRAGMENT INSERTION (SEPP) - Optional for higher accuracy
#-------------------------------------------------------------------------------
# Uncomment below if you want to use fragment insertion with Greengenes/SILVA
# This provides more accurate placement but requires additional reference data

# echo ""
# echo "Running fragment insertion with SEPP (optional, more accurate)..."
# 
# # Download SEPP reference if not present
# SEPP_REF="${DATABASE_DIR}/sepp/sepp-refs-silva-128.qza"
# if [ ! -f "${SEPP_REF}" ]; then
#     echo "Downloading SEPP reference..."
#     wget -O "${SEPP_REF}" \
#         "https://data.qiime2.org/2024.10/common/sepp-refs-silva-128.qza"
# fi
# 
# qiime fragment-insertion sepp \
#     --i-representative-sequences "${INPUT_SEQS}" \
#     --i-reference-database "${SEPP_REF}" \
#     --p-threads ${THREADS} \
#     --o-tree "${OUTPUT_DIR}/sepp-tree.qza" \
#     --o-placements "${OUTPUT_DIR}/sepp-placements.qza"

#-------------------------------------------------------------------------------
# EXPORT TREE
#-------------------------------------------------------------------------------
echo ""
echo "Exporting phylogenetic tree..."

# Export to Newick format
qiime tools export \
    --input-path "${OUTPUT_DIR}/rooted-tree.qza" \
    --output-path "${TABLES_DIR}"

mv "${TABLES_DIR}/tree.nwk" "${TABLES_DIR}/phylogenetic_tree.nwk"

check_status "Tree export"

#-------------------------------------------------------------------------------
# TREE STATISTICS
#-------------------------------------------------------------------------------
echo ""
echo "=============================================="
echo "Phylogenetic Tree Statistics"
echo "=============================================="

# Count tips (should match ASV count)
TIPS=$(grep -o ":" "${TABLES_DIR}/phylogenetic_tree.nwk" | wc -l)
echo "Number of tree tips: $((TIPS/2))"

# Get tree file size
TREE_SIZE=$(ls -lh "${TABLES_DIR}/phylogenetic_tree.nwk" | awk '{print $5}')
echo "Tree file size: ${TREE_SIZE}"

#-------------------------------------------------------------------------------
# SUMMARY
#-------------------------------------------------------------------------------
echo ""
echo "=============================================="
echo "Phylogeny Construction Complete!"
echo "=============================================="
echo "End time: $(date)"
echo ""
echo "Output files:"
echo "  Alignment: ${OUTPUT_DIR}/aligned-rep-seqs.qza"
echo "  Masked alignment: ${OUTPUT_DIR}/masked-aligned-rep-seqs.qza"
echo "  Unrooted tree: ${OUTPUT_DIR}/unrooted-tree.qza"
echo "  Rooted tree: ${OUTPUT_DIR}/rooted-tree.qza"
echo "  Newick export: ${TABLES_DIR}/phylogenetic_tree.nwk"
echo ""
echo "Next step: bash ${SCRIPT_DIR}/09_diversity.sh"

duration=$SECONDS
echo "Phylogeny construction completed in $(($duration / 60)) minutes and $(($duration % 60)) seconds."