#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=80
#SBATCH --time=7-00:00:00
#SBATCH --partition=week
#SBATCH --output=slurm/%x.%j.out

# Script name:  03_extract_subregions.sh
# Description:  Extract ITS1 and ITS2 subregions from the full ITS FASTA using
#               ITSx. Filters the classification table to retain only sequences
#               present in each extracted subregion.
# Note:         This script must be run from the project root directory.

# =============================================================================
# PARAMETER SETUP
# =============================================================================

readonly GLOBAL_FASTA="./data/full_ITS/eukaryome_ITS.fasta"
readonly GLOBAL_CLASS="./data/full_ITS/eukaryome_ITS.classification"

readonly ITSX_TMP="./tmp/itsx"
readonly PREFIX="eukaryome_ITS"

readonly N_CPUS="${SLURM_CPUS_PER_TASK:-$(nproc)}"

# =============================================================================
# DIRECTORY SETUP
# =============================================================================

mkdir -p ./data/ITS1 ./data/ITS2 "$ITSX_TMP"

# =============================================================================
# ENVIRONMENT SETUP
# =============================================================================

echo "Activating conda environment..."
source ~/.bashrc
conda activate dyna_clust_env

# =============================================================================
# INPUT VALIDATION
# =============================================================================

if [[ ! -f "$GLOBAL_FASTA" ]]; then
    echo "ERROR: Full ITS FASTA not found: $GLOBAL_FASTA" >&2
    echo "Please run 01_reformat_ITS.sh first." >&2
    exit 1
fi

if [[ ! -f "$GLOBAL_CLASS" ]]; then
    echo "ERROR: Classification file not found: $GLOBAL_CLASS" >&2
    echo "Please run 01_reformat_ITS.sh first." >&2
    exit 1
fi

# =============================================================================
# RUN ITSx
# =============================================================================

echo ""
echo "=== RUNNING ITSx ==="
echo "$(date)"
echo "Input FASTA : $GLOBAL_FASTA"
echo "Output dir  : $ITSX_TMP"
echo "Threads     : $N_CPUS"
echo ""

ITSx \
    -i  "$GLOBAL_FASTA" \
    -o  "$ITSX_TMP/$PREFIX" \
    --cpu "$N_CPUS" \
    --preserve T \
    -t all \
    --graphical F

if [[ $? -ne 0 ]]; then
    echo "ERROR: ITSx failed." >&2
    exit 1
fi

echo "ITSx completed at: $(date)"
echo ""

# =============================================================================
# HELPER FUNCTION: filter classification to IDs present in a FASTA
# =============================================================================

# filter_classification <fasta> <global_class> <out_class>
filter_classification() {
    local fasta="$1"
    local global_class="$2"
    local out_class="$3"

    # Extract sequence IDs from FASTA headers, stripping the leading '>'
    local id_tmp
    id_tmp=$(mktemp)
    grep "^>" "$fasta" | sed 's/^>//' | awk '{print $1}' > "$id_tmp"

    # Retain header row plus any classification row whose first field matches
    awk 'NR == FNR { ids[$1] = 1; next }
         FNR == 1  { print; next }
         $1 in ids { print }' "$id_tmp" "$global_class" > "$out_class"

    rm -f "$id_tmp"

    local n_rows
    n_rows=$(( $(wc -l < "$out_class") - 1 ))
    echo "  Classification rows written: ${n_rows}"
}

# =============================================================================
# ITS1
# =============================================================================

echo "=== PROCESSING ITS1 ==="

ITS1_RAW="$ITSX_TMP/${PREFIX}.ITS1.fasta"
ITS1_FASTA="./data/ITS1/eukaryome_ITS1.fasta"
ITS1_CLASS="./data/ITS1/eukaryome_ITS1.classification"

if [[ ! -f "$ITS1_RAW" ]]; then
    echo "WARNING: ITSx ITS1 output not found: $ITS1_RAW — skipping ITS1." >&2
else
    cp "$ITS1_RAW" "$ITS1_FASTA"
    ITS1_COUNT=$(grep -c "^>" "$ITS1_FASTA")
    echo "  ITS1 sequences: $ITS1_COUNT"

    echo "  Filtering classification for ITS1..."
    filter_classification "$ITS1_FASTA" "$GLOBAL_CLASS" "$ITS1_CLASS"

    echo "  ITS1 FASTA written to          : $ITS1_FASTA"
    echo "  ITS1 classification written to : $ITS1_CLASS"
fi

echo ""

# =============================================================================
# ITS2
# =============================================================================

echo "=== PROCESSING ITS2 ==="

ITS2_RAW="$ITSX_TMP/${PREFIX}.ITS2.fasta"
ITS2_FASTA="./data/ITS2/eukaryome_ITS2.fasta"
ITS2_CLASS="./data/ITS2/eukaryome_ITS2.classification"

if [[ ! -f "$ITS2_RAW" ]]; then
    echo "WARNING: ITSx ITS2 output not found: $ITS2_RAW — skipping ITS2." >&2
else
    cp "$ITS2_RAW" "$ITS2_FASTA"
    ITS2_COUNT=$(grep -c "^>" "$ITS2_FASTA")
    echo "  ITS2 sequences: $ITS2_COUNT"

    echo "  Filtering classification for ITS2..."
    filter_classification "$ITS2_FASTA" "$GLOBAL_CLASS" "$ITS2_CLASS"

    echo "  ITS2 FASTA written to          : $ITS2_FASTA"
    echo "  ITS2 classification written to : $ITS2_CLASS"
fi

echo ""

# =============================================================================
# SUMMARY
# =============================================================================

echo "=== SEQUENCE COUNT SUMMARY ==="
FULL_COUNT=$(grep -c "^>" "$GLOBAL_FASTA")
echo "  Full ITS sequences : $FULL_COUNT"
[[ -f "$ITS1_FASTA" ]] && echo "  ITS1 sequences     : $(grep -c '^>' "$ITS1_FASTA")"
[[ -f "$ITS2_FASTA" ]] && echo "  ITS2 sequences     : $(grep -c '^>' "$ITS2_FASTA")"
echo ""

# =============================================================================
# CLEANUP
# =============================================================================

echo "=== CLEANUP ==="
echo "Removing ITSx tmp directory: $ITSX_TMP"
rm -rf "$ITSX_TMP"

echo ""
echo "=== PIPELINE COMPLETED SUCCESSFULLY ==="
echo "$(date)"

conda deactivate
