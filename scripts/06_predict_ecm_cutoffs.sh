#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=80
#SBATCH --time=60-00:00:00
#SBATCH --partition=long
#SBATCH --output=logs/%x.%j.out

# Script name:  06b_predict_ecm_cutoffs.sh
# Description:  Predict species-level similarity cutoffs grouped by ECM lineage.
#               Builds ecm_lineage.classification and ecm_lineage.fasta on-the-fly
#               by mapping genera to lineages from ecm_genera.txt and replacing the
#               family column in ecm_family.classification. Rows with no matching
#               genus/lineage pair are dropped. Local species cutoffs are predicted
#               for each lineage with >= 20 sequences and >= 5 subgroups, followed
#               by a global species cutoff prediction.
# Note:         This script must be run from the project root directory.

# =============================================================================
# PARAMETER SETUP
# =============================================================================

readonly PREDICT="./R/predict.R"
readonly TMP_DIR="./tmp"
readonly PREFIX="eukaryome_ecm_lineages"
readonly STEP=0.001
readonly START_THRESH=0.7
readonly END_THRESH=1
readonly MAX_SEQ=100000
readonly MAX_PROP=0.5
readonly MIN_SEQ=20
readonly MIN_GROUPS=5
readonly N_CPUS="${SLURM_CPUS_PER_TASK:-$(nproc)}"
readonly RUN_PARALLEL="yes"

# INPUTS
readonly ECM_GENERA="./data/ecm_genera.txt"
readonly INPUT_FASTA_SRC="./data/full_ITS/ecm_family.fasta"
readonly INPUT_CLASS_SRC="./data/full_ITS/ecm_family.classification"

# OUTPUTS
readonly INPUT_FASTA="./data/full_ITS/ecm_lineage.fasta"
readonly INPUT_CLASS="./data/full_ITS/ecm_lineage.classification"
readonly INPUT_SIM="./data/full_ITS/ecm_lineage.sim"
readonly OUT_DIR="./data/full_ITS"

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

if [[ ! -f "$PREDICT" ]]; then
    echo "ERROR: R script not found: $PREDICT" >&2
    exit 1
fi

for f in "$ECM_GENERA" "$INPUT_FASTA_SRC" "$INPUT_CLASS_SRC"; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: Required input file not found: $f" >&2
        exit 1
    fi
done

# =============================================================================
# BUILD ecm_lineage.classification
# Replace the family column with the lineage name for each genus; drop rows
# whose genus is not present in ecm_genera.txt.
# ecm_genera.txt columns (tab-separated, with header): genus lineage ...
# ecm_family.classification columns (tab-separated):
#   id kingdom phylum class order family genus species
# The family column header is kept as-is so downstream prediction code needs
# no changes. Leading "/" is stripped from lineage values (e.g. "/boletus" -> "boletus").
# =============================================================================

echo ""
echo "=== BUILDING ecm_lineage.classification ==="
echo "$(date)"

awk -F'\t' -v OFS='\t' '
    # Pass 1 (FNR==NR): load genus->lineage map from ecm_genera.txt
    FNR == NR {
        if (FNR == 1) next          # skip header
        gsub(/\r/, "")              # strip Windows CR
        genus   = $1
        lineage = $2
        sub(/^\//, "", lineage)     # strip leading slash
        map[genus] = lineage
        next
    }
    # Pass 2: process ecm_family.classification
    {
        gsub(/\r/, "")              # strip Windows CR
    }
    FNR == 1 {
        print                       # keep header unchanged (family column stays "family")
        next
    }
    {
        genus = $7
        if (genus in map) {
            $6 = map[genus]
            print
        }
        # else: genus not in ECM lineage list — skip row
    }
' "$ECM_GENERA" "$INPUT_CLASS_SRC" > "$INPUT_CLASS"

n_rows=$(awk 'NR > 1' "$INPUT_CLASS" | wc -l)
echo "  Written ${n_rows} sequences to ${INPUT_CLASS}"

# =============================================================================
# BUILD ecm_lineage.fasta
# Keep only FASTA entries whose ID appears in ecm_lineage.classification.
# =============================================================================

echo ""
echo "=== BUILDING ecm_lineage.fasta ==="
echo "$(date)"

awk -F'\t' '
    # Pass 1: collect IDs from classification file (skip header)
    FNR == NR {
        if (FNR == 1) next
        ids[$1] = 1
        next
    }
    # Pass 2: filter FASTA
    /^>/ {
        id = substr($0, 2)
        sub(/[[:space:]].*/, "", id)
        keep = (id in ids) ? 1 : 0
    }
    keep { print }
' "$INPUT_CLASS" "$INPUT_FASTA_SRC" > "$INPUT_FASTA"

n_seqs=$(grep -c "^>" "$INPUT_FASTA")
echo "  Written ${n_seqs} sequences to ${INPUT_FASTA}"

# =============================================================================
# MAIN WORKFLOW
# =============================================================================

echo ""
echo "=== PREDICTING LOCAL SPECIES SIMILARITY CUTOFFS BY LINEAGE ==="
echo "$(date)"

Rscript "$PREDICT" \
    --input           "$INPUT_FASTA" \
    --classification  "$INPUT_CLASS" \
    --sim             "$INPUT_SIM" \
    --rank            species \
    --higher_rank     family \
    --min_seq_no      "$MIN_SEQ" \
    --min_group_no    "$MIN_GROUPS" \
    --max_seq_no      "$MAX_SEQ" \
    --max_proportion  "$MAX_PROP" \
    --start_threshold "$START_THRESH" \
    --end_threshold   "$END_THRESH" \
    --step            "$STEP" \
    --prefix          "$PREFIX" \
    --id_col          id \
    --run_parallel    "$RUN_PARALLEL" \
    --n_cpus          "$N_CPUS" \
    --tmp_dir         "$TMP_DIR" \
    --out             "$OUT_DIR"

echo ""
echo "=== PREDICTING GLOBAL SPECIES SIMILARITY CUTOFFS ==="
echo "$(date)"

Rscript "$PREDICT" \
    --input           "$INPUT_FASTA" \
    --classification  "$INPUT_CLASS" \
    --sim             "$INPUT_SIM" \
    --rank            species \
    --max_seq_no      "$MAX_SEQ" \
    --max_proportion  "$MAX_PROP" \
    --start_threshold "$START_THRESH" \
    --end_threshold   "$END_THRESH" \
    --step            "$STEP" \
    --prefix          "$PREFIX" \
    --id_col          id \
    --run_parallel    "$RUN_PARALLEL" \
    --n_cpus          "$N_CPUS" \
    --tmp_dir         "$TMP_DIR" \
    --out             "$OUT_DIR"

echo ""
echo "Cutoff files written to ${OUT_DIR}:"
echo "  Local cutoffs  : ${PREFIX}.cutoffs.json"
echo "  Global cutoffs : ${PREFIX}.cutoffs.json"

conda deactivate
