#!/bin/bash
# Script name:  06b_launch_all_regions.sh
# Description:  Submit one SLURM job per ITS region (full_ITS, ITS1, ITS2),
#               each with 32 cores, running in parallel instead of sequentially.
#
# Usage:        bash scripts/06b_launch_all_regions.sh
# Note:         Must be run from the project root directory.

set -euo pipefail

SCRIPT="scripts/06b_predict_cutoffs_region.sh"

if [[ ! -f "$SCRIPT" ]]; then
    echo "ERROR: Worker script not found: $SCRIPT" >&2
    exit 1
fi

for region in full_ITS ITS1 ITS2; do
    echo "Submitting ${region}..."
    sbatch --job-name="predict_${region}" "$SCRIPT" "$region"
done

echo ""
echo "All three regions submitted. Monitor with:  squeue -u \$USER"
