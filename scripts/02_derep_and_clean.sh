#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --time=0-01:00:00
#SBATCH --partition=day
#SBATCH --output=slurm/%x.%j.out

# Script name:  02_derep_and_clean.sh
# Description:  Filter, clean, and dereplicate the ITS reference dataset.
#               1. Filter sequences by length (length_filter.R)
#               2. Standardise infraspecific annotations (check_annotations.R)
#               3. Dereplicate identical sequences with LCA taxonomy (dereplicate_lca.R)
# Note:         This script must be run from the project root directory.
# Note:         check_annotations.R needs to be validated for each update of
#               EUKARYOME, as the annotation format may change.

# =============================================================================
# PARAMETER SETUP
# =============================================================================

readonly CLASS="./tmp/eukaryome_ITS.classification"
readonly FASTA="./tmp/eukaryome_ITS.fasta"

readonly CLASS_FILTERED="./tmp/eukaryome_ITS_filtered.classification"
readonly FASTA_FILTERED="./tmp/eukaryome_ITS_filtered.fasta"

readonly CLASS_DEREP="./data/full_ITS/eukaryome_ITS.classification"
readonly FASTA_DEREP="./data/full_ITS/eukaryome_ITS.fasta"

readonly LENGTH_FILTER="./R/length_filter.R"
readonly CHECK_ANNOT="./R/check_annotations.R"
readonly DEREP_LCA="./R/dereplicate_lca.R"

readonly MIN_LENGTH=250
readonly MAX_LENGTH=1500
readonly EXCLUDE_AMBIGUOUS=TRUE  # set to FALSE to keep sequences with ambiguous nucleotides

mkdir -p ./data/full_ITS ./tmp

# =============================================================================
# ENVIRONMENT SETUP
# =============================================================================

echo "Activating conda environment..."
source ~/.bashrc
conda activate dyna_clust_predict

# =============================================================================
# INPUT VALIDATION
# =============================================================================

if [[ ! -f "$CLASS" ]]; then
    echo "ERROR: Classification file not found: $CLASS" >&2
    echo "Please run 01_reformat_ITS.sh first." >&2
    exit 1
fi

if [[ ! -f "$FASTA" ]]; then
    echo "ERROR: FASTA file not found: $FASTA" >&2
    echo "Please run 01_reformat_ITS.sh first." >&2
    exit 1
fi

if [[ ! -f "$LENGTH_FILTER" ]]; then
    echo "ERROR: R script not found: $LENGTH_FILTER" >&2
    exit 1
fi

if [[ ! -f "$CHECK_ANNOT" ]]; then
    echo "ERROR: R script not found: $CHECK_ANNOT" >&2
    exit 1
fi

if [[ ! -f "$DEREP_LCA" ]]; then
    echo "ERROR: R script not found: $DEREP_LCA" >&2
    exit 1
fi

# =============================================================================
# STEP 1: FILTER SEQUENCES BY LENGTH
# =============================================================================

echo ""
echo "=== STEP 1: LENGTH FILTER ==="
echo "$(date)"
echo "FASTA:                         $FASTA"
echo "Classification:                $CLASS"
echo "Min length:                    $MIN_LENGTH"
echo "Max length:                    $MAX_LENGTH"
echo "Exclude ambiguous nucleotides: $EXCLUDE_AMBIGUOUS"
echo ""

AMBIG_FLAG="--exclude_ambiguous"
[[ "$EXCLUDE_AMBIGUOUS" == "FALSE" ]] && AMBIG_FLAG="--no-exclude_ambiguous"

Rscript "$LENGTH_FILTER" \
    --fasta_in             "$FASTA" \
    --fasta_out            "$FASTA_FILTERED" \
    --classification_in    "$CLASS" \
    --classification_out   "$CLASS_FILTERED" \
    --min_length           "$MIN_LENGTH" \
    --max_length           "$MAX_LENGTH" \
    "$AMBIG_FLAG"

if [[ $? -ne 0 ]]; then
    echo "ERROR: length_filter.R failed." >&2
    exit 1
fi

echo ""
echo "Filtered FASTA written to:          $FASTA_FILTERED"
echo "Filtered classification written to: $CLASS_FILTERED"

# =============================================================================
# STEP 2: CHECK AND STANDARDISE INFRASPECIFIC ANNOTATIONS
# =============================================================================

echo ""
echo "=== STEP 2: CHECKING INFRASPECIFIC ANNOTATIONS ==="
echo "$(date)"
echo "Classification file: $CLASS_FILTERED"
echo ""

Rscript "$CHECK_ANNOT" \
    --classification_in  "$CLASS_FILTERED" \
    --classification_out "$CLASS_FILTERED"

if [[ $? -ne 0 ]]; then
    echo "ERROR: check_annotations.R failed." >&2
    exit 1
fi

echo ""
echo "Updated classification written to: $CLASS_FILTERED"

# =============================================================================
# STEP 3: DEREPLICATE AND RESOLVE LCA TAXONOMY
# =============================================================================

echo ""
echo "=== STEP 3: DEREPLICATE AND RESOLVE LCA TAXONOMY ==="
echo "$(date)"
echo "FASTA:           $FASTA_FILTERED"
echo "Classification:  $CLASS_FILTERED"
echo ""

Rscript "$DEREP_LCA" \
    --fasta_in             "$FASTA_FILTERED" \
    --fasta_out            "$FASTA_DEREP" \
    --classification_in    "$CLASS_FILTERED" \
    --classification_out   "$CLASS_DEREP"

if [[ $? -ne 0 ]]; then
    echo "ERROR: dereplicate_lca.R failed." >&2
    exit 1
fi

# Clean the ./tmp directory
rm -rf ./tmp/*

echo ""
echo "Dereplicated FASTA written to:          $FASTA_DEREP"
echo "Dereplicated classification written to: $CLASS_DEREP"
echo ""
echo "=== PIPELINE COMPLETED SUCCESSFULLY ==="
echo "$(date)"

conda deactivate
