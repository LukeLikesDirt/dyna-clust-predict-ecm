#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=3950
#SBATCH --time=60-00:00:00
#SBATCH --partition=long
#SBATCH --nodelist=ltu-hpc-2
#SBATCH --output=slurm/%x.%j.out

# Script name:  06b_predict_cutoffs_region.sh
# Description:  Predict similarity cutoffs for a SINGLE ITS region.
#               Designed to be submitted once per region so full_ITS, ITS1 and
#               ITS2 run as parallel SLURM jobs (each with 32 cores) instead of
#               sequentially in one large job.
#
# Usage:        sbatch --job-name=predict_full_ITS scripts/06b_predict_cutoffs_region.sh full_ITS
#               sbatch --job-name=predict_ITS1     scripts/06b_predict_cutoffs_region.sh ITS1
#               sbatch --job-name=predict_ITS2     scripts/06b_predict_cutoffs_region.sh ITS2
#
#               Or use the launcher:
#               bash scripts/06b_launch_all_regions.sh
#
# Note:         This script must be run from the project root directory.

# =============================================================================
# REGION ARGUMENT
# =============================================================================

region="${1:?ERROR: Region argument required (full_ITS, ITS1, or ITS2)}"

case "$region" in
    full_ITS|ITS1|ITS2) ;;
    *) echo "ERROR: Invalid region '$region'. Must be full_ITS, ITS1, or ITS2." >&2
       exit 1 ;;
esac

# =============================================================================
# PARAMETER SETUP
# =============================================================================

readonly PREDICT="./R/predict.R"
readonly TMP_DIR="./tmp/${region}"      # per-region tmp avoids collisions
readonly PREFIX="eukaryome"
readonly STEP=0.001
readonly END_THRESH=1
readonly N_CPUS="${SLURM_CPUS_PER_TASK:-$(nproc)}"
readonly RUN_PARALLEL="yes"

# Start thresholds per target rank
declare -A START_THRESH
START_THRESH["species"]="0.8"
START_THRESH["genus"]="0.7"
START_THRESH["family"]="0.6"
START_THRESH["order"]="0.5"
START_THRESH["class"]="0.5"
START_THRESH["phylum"]="0.5"

# Parent rank combinations per target rank
declare -A PARENT_RANKS
PARENT_RANKS["species"]="genus family order class phylum kingdom"
PARENT_RANKS["genus"]="family order class phylum kingdom"
PARENT_RANKS["family"]="order class phylum kingdom"
PARENT_RANKS["order"]="class phylum kingdom"
PARENT_RANKS["class"]="phylum kingdom"
PARENT_RANKS["phylum"]="kingdom"

# Filename abbreviations matching subset.R output names
declare -A RANK_ABBR
RANK_ABBR["kingdom"]="kng"
RANK_ABBR["phylum"]="phy"
RANK_ABBR["class"]="cls"
RANK_ABBR["order"]="ord"
RANK_ABBR["family"]="fam"
RANK_ABBR["genus"]="gen"

TARGET_RANKS=("species" "genus" "family" "order" "class" "phylum")

# FASTA and classification paths per region
declare -A REGION_FASTA
REGION_FASTA["full_ITS"]="./data/full_ITS/eukaryome_ITS.fasta"
REGION_FASTA["ITS1"]="./data/ITS1/eukaryome_ITS1.fasta"
REGION_FASTA["ITS2"]="./data/ITS2/eukaryome_ITS2.fasta"

declare -A REGION_CLASS
REGION_CLASS["full_ITS"]="./data/full_ITS/eukaryome_ITS.classification"
REGION_CLASS["ITS1"]="./data/ITS1/eukaryome_ITS1.classification"
REGION_CLASS["ITS2"]="./data/ITS2/eukaryome_ITS2.classification"

# =============================================================================
# DIRECTORY SETUP
# =============================================================================

mkdir -p "$TMP_DIR"

# =============================================================================
# HELPER FUNCTION
# =============================================================================

subset_fasta_and_classification() {
    local id_file="$1"
    local global_fasta="$2"
    local global_class="$3"
    local out_fasta="$4"
    local out_class="$5"

    if [[ ! -f "$id_file" ]]; then
        echo "  WARNING: ID file not found, skipping: $id_file" >&2
        return 1
    fi

    seqkit grep -f "$id_file" "$global_fasta" > "$out_fasta"

    local n_seqs
    n_seqs=$(grep -c "^>" "$out_fasta")
    echo "  Wrote ${n_seqs} sequences to ${out_fasta}"

    awk 'NR == FNR { ids[$1] = 1; next }
         FNR == 1  { print; next }
         $1 in ids { print }' "$id_file" "$global_class" > "$out_class"

    local n_class
    n_class=$(( $(wc -l < "$out_class") - 1 ))
    echo "  Wrote ${n_class} classification rows to ${out_class}"
}

# =============================================================================
# ENVIRONMENT SETUP
# =============================================================================

echo "Activating conda environment..."
source ~/.bashrc
conda activate dyna_clust_predict

# Track failures for the end-of-job summary
declare -a FAILED_STEPS=()

# =============================================================================
# INPUT VALIDATION
# =============================================================================

if [[ ! -f "$PREDICT" ]]; then
    echo "ERROR: R script not found: $PREDICT" >&2
    exit 1
fi

global_fasta="${REGION_FASTA[$region]}"
global_class="${REGION_CLASS[$region]}"

if [[ ! -f "$global_fasta" ]]; then
    echo "ERROR: FASTA not found for region '$region': $global_fasta" >&2
    exit 1
fi

if [[ ! -f "$global_class" ]]; then
    echo "ERROR: Classification not found for region '$region': $global_class" >&2
    exit 1
fi

# =============================================================================
# MAIN: process a single region
# =============================================================================

echo ""
echo "=== REGION: $region ==="
echo "$(date)"

# -------------------------------------------------------------------------
# LOCAL PREDICTIONS (per target x parent combination)
# -------------------------------------------------------------------------

echo ""
echo "--- LOCAL PREDICTIONS: $region ---"
echo "$(date)"

for target in "${TARGET_RANKS[@]}"; do

    st="${START_THRESH[$target]}"

    for parent in ${PARENT_RANKS[$target]}; do

        parent_abbr="${RANK_ABBR[$parent]}"
        id_file="./data/${region}/${target}_pred_id_${parent_abbr}.txt"
        tmp_fasta="$TMP_DIR/${target}_pred_${parent_abbr}.fasta"
        tmp_class="$TMP_DIR/${target}_pred_${parent_abbr}.classification"

        if [[ ! -f "$id_file" ]]; then
            echo "  Skipping ${target} within ${parent}: ID file not found ($id_file)"
            continue
        fi

        echo ""
        echo "  -- ${region} / ${target} within ${parent} --"
        echo "  Subsetting sequences at: $(date)"
        subset_fasta_and_classification \
            "$id_file" "$global_fasta" "$global_class" \
            "$tmp_fasta" "$tmp_class" || continue

        echo "  Predicting cutoffs at: $(date)"
        Rscript "$PREDICT" \
            --input          "$tmp_fasta" \
            --classification "$tmp_class" \
            --rank           "$target" \
            --higher_rank    "$parent" \
            --start_threshold "${st}" \
            --end_threshold  "$END_THRESH" \
            --step           "$STEP" \
            --prefix         "$PREFIX" \
            --id_col         id \
            --run_parallel   "$RUN_PARALLEL" \
            --n_cpus         "$N_CPUS" \
            --tmp_dir        "$TMP_DIR" \
            --out            "./data/${region}"

        rc=$?
        if [[ $rc -ne 0 ]]; then
            echo "" >&2
            echo "  ============================================================" >&2
            echo "  FAILED: ${region}/${target} within ${parent}  (exit code $rc)" >&2
            echo "  Time of failure: $(date)" >&2
            echo "  Memory at failure (RSS peak of this job):" >&2
            echo "    $(grep VmPeak /proc/$$/status 2>/dev/null || echo 'N/A')" >&2
            echo "    $(grep VmRSS  /proc/$$/status 2>/dev/null || echo 'N/A')" >&2
            if command -v sstat &>/dev/null && [[ -n "${SLURM_JOB_ID:-}" ]]; then
                echo "  SLURM memory usage:" >&2
                sstat -j "${SLURM_JOB_ID}.batch" --format=JobID,MaxRSS,MaxVMSize 2>/dev/null || true
            fi
            echo "  If the R output above mentions 'FutureInterruptError' or" >&2
            echo "  'parallel job did not deliver a result', the likely cause" >&2
            echo "  is the OS OOM killer terminating worker processes." >&2
            echo "  ============================================================" >&2
            echo "" >&2
            FAILED_STEPS+=("${target} within ${parent}")
        else
            echo "  Finished at: $(date)"
        fi

        rm -f "$tmp_fasta" "$tmp_class"

    done
done

# -------------------------------------------------------------------------
# GLOBAL PREDICTIONS (no higher rank, per target rank)
# -------------------------------------------------------------------------

echo ""
echo "--- GLOBAL PREDICTIONS: $region ---"
echo "$(date)"

for target in "${TARGET_RANKS[@]}"; do

    st="${START_THRESH[$target]}"
    id_file="./data/${region}/${target}_unique_id.txt"
    tmp_fasta="$TMP_DIR/${target}_unique.fasta"
    tmp_class="$TMP_DIR/${target}_unique.classification"

    if [[ ! -f "$id_file" ]]; then
        echo "  Skipping global prediction for ${target}: ID file not found ($id_file)"
        continue
    fi

    echo ""
    echo "  -- ${region} / global: ${target} --"
    echo "  Subsetting sequences at: $(date)"
    subset_fasta_and_classification \
        "$id_file" "$global_fasta" "$global_class" \
        "$tmp_fasta" "$tmp_class" || continue

    echo "  Predicting global cutoffs at: $(date)"
    Rscript "$PREDICT" \
        --input          "$tmp_fasta" \
        --classification "$tmp_class" \
        --rank           "$target" \
        --start_threshold "${st}" \
        --end_threshold  "$END_THRESH" \
        --step           "$STEP" \
        --prefix         "$PREFIX" \
        --id_col         id \
        --run_parallel   "$RUN_PARALLEL" \
        --n_cpus         "$N_CPUS" \
        --tmp_dir        "$TMP_DIR" \
        --out            "./data/${region}"

    rc=$?
    if [[ $rc -ne 0 ]]; then
        echo "" >&2
        echo "  ============================================================" >&2
        echo "  FAILED: ${region}/global ${target}  (exit code $rc)" >&2
        echo "  Time of failure: $(date)" >&2
        echo "  ============================================================" >&2
        echo "" >&2
        FAILED_STEPS+=("global ${target}")
    else
        echo "  Finished at: $(date)"
    fi

    rm -f "$tmp_fasta" "$tmp_class"

done

# =============================================================================
# CLEANUP
# =============================================================================

echo ""
echo "=== CLEANUP ==="
echo "Removing tmp directory: $TMP_DIR"
rm -rf "$TMP_DIR"

echo ""
if [[ ${#FAILED_STEPS[@]} -gt 0 ]]; then
    echo "=== REGION $region COMPLETED WITH ${#FAILED_STEPS[@]} FAILURE(S) ==="
    echo "$(date)"
    echo ""
    echo "Failed steps:"
    for step in "${FAILED_STEPS[@]}"; do
        echo "  - ${region} / ${step}"
    done
    echo ""
    echo "Likely cause: OS OOM (out-of-memory) killer terminated worker processes."
    echo "Possible fixes:"
    echo "  1. Request more memory in SLURM (#SBATCH --mem=...)."
    echo "  2. Reduce --max_seq_no (currently ${max_seq_no:-25000}) to lower per-dataset memory."
    echo "  3. Reduce --n_cpus to run fewer parallel workers."
else
    echo "=== REGION $region COMPLETED SUCCESSFULLY ==="
    echo "$(date)"
fi
echo ""
echo "Cutoff files written to data/${region}"
echo "  Local cutoffs  : ${PREFIX}.cutoffs.json"
echo "  Global cutoffs : ${PREFIX}.cutoffs.json"

conda deactivate
