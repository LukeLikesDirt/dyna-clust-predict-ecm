#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH --mem-per-cpu=4G
#SBATCH --time=1-00:00:00
#SBATCH --partition=day
#SBATCH --output=logs/%x.%j.out

# Script name:  05b_compute_species_sim.sh
# Description:  Pre-compute a pairwise vsearch similarity matrix for all
#               sequences identified at species rank, for each ITS region.
#               Species-identified IDs are extracted directly from the
#               classification file (matching the is_identified logic in utils.R)
#               and used to subset the full FASTA before running compute_sim.R.
#
#               Outputs (one per region):
#                 ./data/full_ITS/species_ITS.sim
#                 ./data/ITS1/species_ITS1.sim
#                 ./data/ITS2/species_ITS2.sim
#
#               Skips any region whose .sim file already exists.
# Note:         This script must be run from the project root directory.

# =============================================================================
# PARAMETER SETUP
# =============================================================================

readonly SIM="./R/compute_sim.R"
readonly TMP_DIR="./tmp"
readonly MIN_SIM=0.7
readonly N_CPUS="${SLURM_CPUS_PER_TASK:-$(nproc)}"
readonly INPUT_FASTA="data/full_ITS/ecm_family.fasta"
readonly OUTPUT_DIR="data/full_ITS/"

# =============================================================================
# DIRECTORY SETUP
# =============================================================================

mkdir -p "$TMP_DIR"

# =============================================================================
# ENVIRONMENT SETUP
# =============================================================================

echo "Activating conda environment..."
source ~/.bashrc
conda activate dyna_clust_predict

# =============================================================================
# INPUT VALIDATION
# =============================================================================

if [[ ! -f "$SIM" ]]; then
    echo "ERROR: R script not found: $SIM" >&2
    exit 1
fi

# =============================================================================
# COMPUTE SIMILARITY MATRIx
# ============================================================================

echo ""
echo "=== COMPUTING VSEARCH SPECIES SIMILARITY MATRICES ==="
echo "$(date)"
echo "Threads : $N_CPUS"
echo "min_sim : $MIN_SIM"
echo ""
Rscript "$SIM" \
        --input    "$INPUT_FASTA" \
        --out      "$OUTPUT_DIR" \
        --min_sim  "$MIN_SIM" \
        --n_cpus   "$N_CPUS"

# =============================================================================
# CLEANUP
# =============================================================================

echo "=== CLEANUP ==="
echo "Removing tmp directory..."
rm -rf "$TMP_DIR"

echo ""
echo "=== PIPELINE COMPLETED SUCCESSFULLY ==="
echo "$(date)"
echo ""
echo "Sim files written:"


conda deactivate
