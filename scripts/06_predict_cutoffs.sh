#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=80
#SBATCH --time=60-00:00:00
#SBATCH --partition=long
#SBATCH --output=slurm/%x.%j.out

# Script name:  06_predict_cutoffs.sh
# Description:  Predict similarity cutoffs for all three ITS regions (full_ITS,
#               ITS1, ITS2) sequentially. Within each region, local predictions
#               are run for every (target x parent) rank combination followed by
#               global predictions (no higher rank). Sequence subset FASTA and
#               classification files are built on-the-fly from each region's
#               master FASTA to keep tmp/ lean.
# Note:         This script must be run from the project root directory.

# =============================================================================
# PARAMETER SETUP
# =============================================================================

readonly PREDICT="./R/predict.R"
readonly TMP_DIR="./tmp"
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
REGIONS=("full_ITS" "ITS1" "ITS2")

# FASTA path per region
declare -A REGION_FASTA
REGION_FASTA["full_ITS"]="./data/full_ITS/eukaryome_ITS.fasta"
REGION_FASTA["ITS1"]="./data/ITS1/eukaryome_ITS1.fasta"
REGION_FASTA["ITS2"]="./data/ITS2/eukaryome_ITS2.fasta"

# Classification path per region
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

# subset_fasta_and_classification <id_file> <global_fasta> <global_class>
#                                 <out_fasta> <out_class>
# Subsets a global FASTA and classification table to only the IDs listed in
# id_file. Uses seqkit for FASTA filtering and a two-pass awk for the TSV.
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

    # Filter FASTA by ID list using seqkit grep.
    # -f reads IDs from a file (one per line); output preserves full headers.
    seqkit grep -f "$id_file" "$global_fasta" > "$out_fasta"

    local n_seqs
    n_seqs=$(grep -c "^>" "$out_fasta")
    echo "  Wrote ${n_seqs} sequences to ${out_fasta}"

    # Filter classification TSV: retain the header row plus any row whose
    # first field (sequence ID) appears in the ID file (two-pass awk).
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

# =============================================================================
# INPUT VALIDATION
# =============================================================================

if [[ ! -f "$PREDICT" ]]; then
    echo "ERROR: R script not found: $PREDICT" >&2
    exit 1
fi

# =============================================================================
# MAIN LOOP: iterate over all three regions
# =============================================================================

for region in "${REGIONS[@]}"; do

    global_fasta="${REGION_FASTA[$region]}"
    global_class="${REGION_CLASS[$region]}"

    echo ""
    echo "=== REGION: $region ==="
    echo "$(date)"

    if [[ ! -f "$global_fasta" ]]; then
        echo "WARNING: FASTA not found, skipping region '$region': $global_fasta" >&2
        continue
    fi

    if [[ ! -f "$global_class" ]]; then
        echo "WARNING: Classification not found, skipping region '$region': $global_class" >&2
        continue
    fi

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

            # Skip if the ID file was not produced (no groups passed filters)
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

            if [[ $? -ne 0 ]]; then
                echo "  WARNING: predict.R failed for ${region}/${target} within ${parent}" >&2
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

        if [[ $? -ne 0 ]]; then
            echo "  WARNING: global predict.R failed for ${region}/${target}" >&2
        else
            echo "  Finished at: $(date)"
        fi

        rm -f "$tmp_fasta" "$tmp_class"

    done

    echo ""
    echo "=== REGION $region COMPLETE at $(date) ==="

done

# =============================================================================
# CLEANUP
# =============================================================================

echo ""
echo "=== CLEANUP ==="
echo "Removing tmp directory..."
rm -rf "$TMP_DIR"

echo ""
echo "=== PIPELINE COMPLETED SUCCESSFULLY ==="
echo "$(date)"
echo ""
echo "Cutoff files written to data/full_ITS, data/ITS1, data/ITS2"
echo "  Local cutoffs  : ${PREFIX}.cutoffs.json (per region)"
echo "  Global cutoffs : ${PREFIX}.cutoffs.json (appended per region)"

conda deactivate
