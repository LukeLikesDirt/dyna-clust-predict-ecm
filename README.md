# dyna-clust-predict

Prediction of optimal sequence similarity cut-offs for classification and clustering of metabarcoding data using vsearch global alignment and F-measure optimisation as a confidence measure, following Vu *et al.* (2022).

## Overview

The pipeline takes a curated ITS reference database (EUKARYOME), extracts ITS1 and ITS2 subregions with ITSx, and predicts optimal similarity thresholds at each taxonomic rank for three sequence regions: full ITS, ITS1, and ITS2. This pipeline is conceptually similar to and adapted from [dnabarcoder](https://github.com/vuthuyduong/dnabarcoder), but uses vsearch for global alignment and supports parallel processing to speed up similarity computations and cut-off predictions within an R environment.

## Pipeline steps

| Script | Description |
|--------|-------------|
| `scripts/01_reformat_ITS.sh` | Download EUKARYOME, reformat headers, extract taxonomy |
| `scripts/02_check_annotations.sh` | Standardise infraspecific annotations in classification |
| `scripts/03_extract_subregions.sh` | Run ITSx to extract ITS1 and ITS2 subregions |
| `scripts/04_prepare_subsets.sh` | Generate prediction ID files for each region |
| `scripts/05_compute_sim.sh` | **Optional**: Pre-compute similarity matrices |
| `scripts/06_predict_cut-offs.sh` | Predict cut-offs for full ITS, ITS1, and ITS2 |

All scripts must be run from the **project root directory**.

## Data directory

```
data/
├── full_ITS/    # Full ITS sequences, classification, ID files, predictions
├── ITS1/        # ITS1 subregion sequences, classification, ID files, predictions
└── ITS2/        # ITS2 subregion sequences, classification, ID files, predictions
```

## R scripts

| Script | Description |
|--------|-------------|
| `R/utils.R` | Shared utility: `is_identified()` for taxonomy filtering |
| `R/reformat.R` | Reformat raw EUKARYOME FASTA headers and extract taxonomy |
| `R/check_annotations.R` | Detect and standardise infraspecific annotations |
| `R/subset.R` | Generate balanced prediction subsets and ID files |
| `R/compute_sim.R` | Pre-compute pairwise vsearch similarity matrix |
| `R/predict.R` | Predict optimal similarity cut-offs (parallel or sequential) |

## Environment setup

A [conda](https://docs.conda.io/) environment specification is provided in `environment.yml`.
[mamba](https://mamba.readthedocs.io/) is strongly recommended over `conda` for environment creation due to faster dependency resolution

**Local setup:**
```bash
mamba env create -f environment.yml
conda activate dyna_clust_predict
```

**HPC (SLURM) setup** — run from the project root:
```bash
sbatch scripts/create_dyna_clust_env.sh
conda activate dyna_clust_predict
```

## Quick start

```bash
# Activate environment, then from the project root:
conda activate dyna_clust_predict

sbatch scripts/01_reformat_ITS.sh
sbatch scripts/02_check_annotations.sh
sbatch scripts/03_extract_subregions.sh
sbatch scripts/04_prepare_subsets.sh
sbatch scripts/05_compute_sim.sh   # optional; similarity is computed on the fly in step 06, which is the preferred for larger datasets and parallel processing
sbatch scripts/06_predict_cut-offs.sh
```

## Key parameters for sequence selection: `subset.R` run via `04_prepare_subsets.sh`

| Argument | Default | Description |
|----------|---------|-------------|
| `--min_subgroups` | 10 | Minimum number of unique child taxa required per parent taxon |
| `--min_sequences` | 30 | Minimum number of sequences per parent taxon after proportion cap |
| `--max_sequences` | 25000 | Maximum number of sequences per parent taxon; excess is downsampled proportionally across child taxa |
| `--max_proportion` | 0.5 | Maximum fraction a single child taxon may represent within a parent taxon |

## Key parameters for similarity prediction: `predict.R` run via `06_predict_cut-offs.sh`

| Argument | Default | Description |
|----------|---------|-------------|
| `--start_threshold` | 0 | Starting similarity threshold |
| `--end_threshold` | 1 | Ending similarity threshold |
| `--step` | 0.001 | Threshold step size |
| `--run_parallel` | yes | Enable parallel dataset processing (yes/no) |
| `--n_cpus` | all-1 | Number of CPUs available for parallel workers and vsearch threads |
| `--tmp_dir` | ./tmp | Directory for temporary vsearch output |

## Citation
Duong Vu, R. Henrik Nilsson, Gerard J.M. Verkley (2022). dnabarcoder: an open-source software package for analyzing and predicting DNA sequence similarity cut-offs for fungal sequence identification. Molecular Ecology Resources. https://doi.org/10.1111/1755-0998.13651