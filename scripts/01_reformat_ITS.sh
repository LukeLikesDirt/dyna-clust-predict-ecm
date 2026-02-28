#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --time=1-00:00:00
#SBATCH --partition=day
#SBATCH --output=slurm/%x.%j.out

# Script name:  01_reformat_ITS.sh
# Description:  Download the EUKARYOME ITS v2.0 database, reformat sequence
#               headers with reformat.R, and extract the taxonomy classification
#               table.
# Note:         This script must be run from the project root directory.

# =============================================================================
# PARAMETER SETUP
# =============================================================================

# EUKARYOME download URL (EUKARYOME General EUK ITS v2.0)
readonly EUKARYOME_URL="https://sisu.ut.ee/wp-content/uploads/sites/643/General_EUK_ITS_v2.0.zip"
readonly DOWNLOAD_FILE="./tmp/General_EUK_ITS_v2.0.zip"
readonly EXTRACTED_FASTA="./tmp/General_EUK_ITS_v2.0.fasta"

# Output files
readonly OUT_FASTA="./data/full_ITS/eukaryome_ITS.fasta"
readonly OUT_CLASS="./data/full_ITS/eukaryome_ITS.classification"

# Helper R script
readonly REFORMAT="./R/reformat.R"

# =============================================================================
# DIRECTORY SETUP
# =============================================================================

mkdir -p ./data/full_ITS tmp

# =============================================================================
# ENVIRONMENT SETUP
# =============================================================================

echo "Activating conda environment..."
source ~/.bashrc
conda activate dyna_clust_predict

# =============================================================================
# INPUT VALIDATION
# =============================================================================

if [[ ! -f "$REFORMAT" ]]; then
    echo "ERROR: R script not found: $REFORMAT" >&2
    exit 1
fi

# =============================================================================
# FILE DOWNLOAD: EUKARYOME
# =============================================================================

echo ""
echo "=== DOWNLOADING EUKARYOME DATABASE ==="
echo "$(date)"
echo "URL: $EUKARYOME_URL"

if ! curl -L -o "$DOWNLOAD_FILE" "$EUKARYOME_URL"; then
    echo "ERROR: Failed to download file from $EUKARYOME_URL" >&2
    exit 1
fi

echo "Unzipping downloaded file..."
if ! 7z x "$DOWNLOAD_FILE" -o"./tmp/" -y; then
    echo "ERROR: Failed to unzip $DOWNLOAD_FILE" >&2
    exit 1
fi

echo "Download and extraction completed successfully."
echo ""

# =============================================================================
# REFORMAT HEADERS
# =============================================================================

echo "=== REFORMATTING HEADERS ==="
echo "$(date)"

if [[ ! -f "$EXTRACTED_FASTA" ]]; then
    echo "ERROR: Extracted FASTA not found: $EXTRACTED_FASTA" >&2
    exit 1
fi

echo "Running reformat.R..."
Rscript "$REFORMAT" \
    --fasta_in        "$EXTRACTED_FASTA" \
    --fasta_out       "$OUT_FASTA" \
    --classification_out "$OUT_CLASS"

if [[ $? -ne 0 ]]; then
    echo "ERROR: reformat.R failed." >&2
    exit 1
fi

echo ""
echo "Reformatted FASTA written to      : $OUT_FASTA"
echo "Classification table written to   : $OUT_CLASS"
echo ""

# =============================================================================
# CLEANUP
# =============================================================================

echo "=== CLEANUP ==="
echo "Removing downloaded and extracted files from tmp/..."
rm -f "$DOWNLOAD_FILE" "$EXTRACTED_FASTA"

echo ""
echo "=== PIPELINE COMPLETED SUCCESSFULLY ==="
echo "$(date)"

conda deactivate
