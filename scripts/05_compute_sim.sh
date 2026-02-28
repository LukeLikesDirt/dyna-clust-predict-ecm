#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=80
#SBATCH --time=60-00:00:00
#SBATCH --partition=long
#SBATCH --output=slurm/%x.%j.out

# Script name:  05_compute_sim.sh
# Description:  Optional pre-computation of pairwise vsearch similarity matrices
#               for global predictions. Runs compute_sim.R for each rank's
#               unique-sequence FASTA across all three ITS regions (full_ITS,
#               ITS1, ITS2). Skips any rank whose .sim file already exists.
#               The resulting .sim files can be passed to predict.R via --sim
#               to avoid re-running vsearch during prediction.
#
#               Note: sim computation is O(N²) in sequence pairs. Lower ranks
#               (species, genus) can take many hours; higher ranks complete
#               faster. Comment out entries in TARGET_RANKS if not needed.
# Note:         This script must be run from the project root directory.

# =============================================================================
# PARAMETER SETUP
# =============================================================================

readonly SIM="./R/compute_sim.R"
readonly TMP_DIR="./tmp"
readonly MIN_SIM=0
readonly N_CPUS="${SLURM_CPUS_PER_TASK:-$(nproc)}"

TARGET_RANKS=("species" "genus" "family" "order" "class" "phylum")
REGIONS=("full_ITS" "ITS1" "ITS2")

# =============================================================================
# DIRECTORY SETUP
# =============================================================================

mkdir -p "$TMP_DIR" slurm

# =============================================================================
# HELPER FUNCTION: resolve FASTA path for a given region
# =============================================================================

region_fasta() {
    local region="$1"
    case "$region" in
        full_ITS) echo "./data/full_ITS/eukaryome_ITS.fasta"  ;;
        ITS1)     echo "./data/ITS1/eukaryome_ITS1.fasta"     ;;
        ITS2)     echo "./data/ITS2/eukaryome_ITS2.fasta"     ;;
        *)
            echo "ERROR: Unknown region: $region" >&2
            return 1
            ;;
    esac
}

# =============================================================================
# ENVIRONMENT SETUP
# =============================================================================

echo "Activating conda environment..."
source ~/.bashrc
conda activate dyna_clust_env

# =============================================================================
# INPUT VALIDATION
# =============================================================================

if [[ ! -f "$SIM" ]]; then
    echo "ERROR: R script not found: $SIM" >&2
    exit 1
fi

# =============================================================================
# COMPUTE SIMILARITY MATRICES
# =============================================================================

echo ""
echo "=== COMPUTING VSEARCH SIMILARITY MATRICES ==="
echo "$(date)"
echo "Threads : $N_CPUS"
echo "min_sim : $MIN_SIM"
echo ""

for region in "${REGIONS[@]}"; do

    global_fasta=$(region_fasta "$region")

    echo "=== REGION: $region ==="

    if [[ ! -f "$global_fasta" ]]; then
        echo "  WARNING: FASTA not found, skipping region '$region': $global_fasta" >&2
        echo ""
        continue
    fi

    for target in "${TARGET_RANKS[@]}"; do

        id_file="./data/${region}/${target}_unique_id.txt"
        sim_out="./data/${region}/${target}_unique.sim"
        tmp_fasta="$TMP_DIR/${region}_${target}_unique.fasta"

        echo "--- ${region} / ${target} ---"

        # Skip if sim file already exists
        if [[ -f "$sim_out" ]]; then
            echo "  Sim file already exists, skipping: $sim_out"
            continue
        fi

        # Skip if ID file is missing
        if [[ ! -f "$id_file" ]]; then
            echo "  Skipping: ID file not found ($id_file)"
            continue
        fi

        echo "  Subsetting sequences at: $(date)"
        seqkit grep -f "$id_file" "$global_fasta" > "$tmp_fasta"

        n_seqs=$(grep -c "^>" "$tmp_fasta")
        echo "  Sequences in subset: $n_seqs"

        echo "  Computing similarity at: $(date)"
        Rscript "$SIM" \
            --input  "$tmp_fasta" \
            --out    "./data/${region}" \
            --minsim "$MIN_SIM" \
            --ncpus  "$N_CPUS"

        if [[ $? -ne 0 ]]; then
            echo "  WARNING: compute_sim.R failed for ${region}/${target}" >&2
        else
            # Rename output to the rank-level convention if compute_sim.R named
            # the file after the tmp FASTA basename rather than the sim_out name.
            local_sim="./data/${region}/${region}_${target}_unique.sim"
            if [[ -f "$local_sim" && "$local_sim" != "$sim_out" ]]; then
                mv "$local_sim" "$sim_out"
            fi
            echo "  Finished at: $(date)"
            echo "  Sim file: $sim_out"
        fi

        rm -f "$tmp_fasta"
        echo ""

    done

    echo ""

done

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
echo "Sim files written per region:"
for region in "${REGIONS[@]}"; do
    for target in "${TARGET_RANKS[@]}"; do
        sim_out="./data/${region}/${target}_unique.sim"
        [[ -f "$sim_out" ]] && echo "  $sim_out"
    done
done

conda deactivate
