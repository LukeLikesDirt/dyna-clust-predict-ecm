#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --time=1-00:00:00
#SBATCH --partition=day
#SBATCH --output=slurm/%x.%j.out
#
# Create the dyna_clust_predict conda environment.
# Must be submitted from the project root directory:
#   sbatch env/create_dyna_clust_env.sh

ENV_FILE="environment.yml"
ENV_NAME="dyna_clust_predict"

if ! conda info --envs | grep -q "${ENV_NAME}"; then
  echo "Creating conda environment '${ENV_NAME}' at: $(date)"
  mamba env create -f "${ENV_FILE}"
else
  echo "Environment '${ENV_NAME}' already exists."
fi

