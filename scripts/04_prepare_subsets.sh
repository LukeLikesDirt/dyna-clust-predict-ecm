#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --time=0-02:00:00
#SBATCH --partition=day
#SBATCH --output=slurm/%x.%j.out

# Script name:  04_prepare_subsets.sh
# Description:  Prepare prediction ID files for all three ITS regions (full ITS,
#               ITS1, ITS2) by running subset.R once per region. Produces one
#               ID file per unique-sequence rank (STEP 1) and one ID file per
#               valid (target x parent) rank combination (STEP 2).
# Note:         This script must be run from the project root directory.

# =============================================================================
# PARAMETER SETUP
# =============================================================================

readonly SUBSET="./R/subset.R"

# Filter constants
readonly MIN_SUBGROUPS=10
readonly MIN_SEQUENCES=30
readonly MAX_SEQUENCES=25000
readonly MAX_PROPORTION=0.5

# =============================================================================
# ENVIRONMENT SETUP
# =============================================================================

echo "Activating conda environment..."
source ~/.bashrc
conda activate dyna_clust_predict

# =============================================================================
# INPUT VALIDATION
# =============================================================================

if [[ ! -f "$SUBSET" ]]; then
    echo "ERROR: R script not found: $SUBSET" >&2
    exit 1
fi

# =============================================================================
# HELPER FUNCTION
# =============================================================================

# run_subset <fasta_in> <classification_in> <output_dir> <region_label>
run_subset() {
    local fasta_in="$1"
    local classification_in="$2"
    local output_dir="$3"
    local label="$4"

    if [[ ! -f "$fasta_in" ]]; then
        echo "WARNING: FASTA not found for region '$label', skipping: $fasta_in" >&2
        return 0
    fi

    if [[ ! -f "$classification_in" ]]; then
        echo "WARNING: Classification not found for region '$label', skipping: $classification_in" >&2
        return 0
    fi

    echo ""
    echo "--- Region: $label ---"
    echo "FASTA          : $fasta_in"
    echo "Classification : $classification_in"
    echo "Output dir     : $output_dir"
    echo "$(date)"

    Rscript "$SUBSET" \
        --fasta_in          "$fasta_in" \
        --classification_in "$classification_in" \
        --output_dir        "$output_dir" \
        --min_subgroups     "$MIN_SUBGROUPS" \
        --min_sequences     "$MIN_SEQUENCES" \
        --max_sequences     "$MAX_SEQUENCES" \
        --max_proportion    "$MAX_PROPORTION"

    if [[ $? -ne 0 ]]; then
        echo "ERROR: subset.R failed for region '$label'." >&2
        return 1
    fi

    echo "Finished region '$label' at: $(date)"
}

# =============================================================================
# PREPARE SUBSETS FOR ALL THREE REGIONS
# =============================================================================

echo ""
echo "=== PREPARING SUBSETS ==="
echo "min_subgroups  : $MIN_SUBGROUPS"
echo "min_sequences  : $MIN_SEQUENCES"
echo "max_sequences  : $MAX_SEQUENCES"
echo "max_proportion : $MAX_PROPORTION"

# 1. Full ITS
run_subset \
    "./data/full_ITS/eukaryome_ITS.fasta" \
    "./data/full_ITS/eukaryome_ITS.classification" \
    "./data/full_ITS" \
    "full_ITS"

# 2. ITS1
run_subset \
    "./data/ITS1/eukaryome_ITS1.fasta" \
    "./data/ITS1/eukaryome_ITS1.classification" \
    "./data/ITS1" \
    "ITS1"

# 3. ITS2
run_subset \
    "./data/ITS2/eukaryome_ITS2.fasta" \
    "./data/ITS2/eukaryome_ITS2.classification" \
    "./data/ITS2" \
    "ITS2"

echo ""
echo "=== PIPELINE COMPLETED SUCCESSFULLY ==="
echo "$(date)"
echo ""
echo "ID files written to data/full_ITS, data/ITS1, data/ITS2"
echo "  STEP 1 unique-sequence IDs : <rank>_unique_id.txt"
echo "  STEP 2 prediction IDs      : <target>_pred_id_<parent>.txt"

conda deactivate
