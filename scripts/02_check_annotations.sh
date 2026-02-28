#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --time=0-01:00:00
#SBATCH --partition=day
#SBATCH --output=slurm/%x.%j.out

# Script name:  02_check_annotations.sh
# Description:  Standardise infraspecific annotations in the classification
#               file produced by 01_reformat_ITS.sh. The R script reads and
#               overwrites the same classification file in-place.
# Note:         This script must be run from the project root directory.
# Note:         check_annotations.R needs to be validated for each update of 
#               EUKARYOME, as the annotation format may change.

# =============================================================================
# PARAMETER SETUP
# =============================================================================

readonly CLASS="./data/full_ITS/eukaryome_ITS.classification"
readonly CHECK_ANNOT="./R/check_annotations.R"

# =============================================================================
# ENVIRONMENT SETUP
# =============================================================================

echo "Activating conda environment..."
source ~/.bashrc
conda activate dyna_clust_env

# =============================================================================
# INPUT VALIDATION
# =============================================================================

if [[ ! -f "$CLASS" ]]; then
    echo "ERROR: Classification file not found: $CLASS" >&2
    echo "Please run 01_reformat_ITS.sh first." >&2
    exit 1
fi

if [[ ! -f "$CHECK_ANNOT" ]]; then
    echo "ERROR: R script not found: $CHECK_ANNOT" >&2
    exit 1
fi

# =============================================================================
# CHECK AND STANDARDISE INFRASPECIFIC ANNOTATIONS
# =============================================================================

echo ""
echo "=== CHECKING INFRASPECIFIC ANNOTATIONS ==="
echo "$(date)"
echo "Classification file: $CLASS"
echo ""

Rscript "$CHECK_ANNOT" \
    --classification_in  "$CLASS" \
    --classification_out "$CLASS"

if [[ $? -ne 0 ]]; then
    echo "ERROR: check_annotations.R failed." >&2
    exit 1
fi

echo ""
echo "Updated classification written to: $CLASS"
echo ""
echo "=== PIPELINE COMPLETED SUCCESSFULLY ==="
echo "$(date)"

conda deactivate
