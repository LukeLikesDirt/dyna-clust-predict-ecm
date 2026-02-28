# dyna-clust-predict

Prediction of optimal sequence similarity cut-offs for classification
and clustering of metabarcoding data using **vsearch global alignment**
and **F-measure optimisation** as a confidence metric, following Vu *et
al.* (2022).

------------------------------------------------------------------------

## Overview

`dyna-clust-predict` predicts optimal sequence similarity thresholds
across taxonomic ranks using a curated ITS reference database
(EUKARYOME).

The pipeline:

-   Extracts ITS1 and ITS2 subregions using ITSx\
-   Computes global pairwise similarity using vsearch\
-   Optimises similarity cut-offs using F-measure\
-   Supports parallel processing for efficient large-scale analyses

The workflow is conceptually adapted from
[dnabarcoder](https://github.com/vuthuyduong/dnabarcoder), but:

-   Uses **vsearch** for global alignment\
-   Supports scalable parallel computation\
-   Integrates similarity prediction within an R-based workflow

------------------------------------------------------------------------

## Pipeline structure

All scripts must be run from the **project root directory**.

### Shell scripts

  --------------------------------------------------------------------------------
  Script                               Description
  ------------------------------------ -------------------------------------------
  `scripts/01_reformat_ITS.sh`         Download EUKARYOME, reformat headers,
                                       extract taxonomy

  `scripts/02_check_annotations.sh`    Standardise infraspecific annotations

  `scripts/03_extract_subregions.sh`   Extract ITS1 and ITS2 using ITSx

  `scripts/04_prepare_subsets.sh`      Generate balanced prediction subsets

  `scripts/05_compute_sim.sh`          *(Optional)* Pre-compute similarity
                                       matrices

  `scripts/06_predict_cut-offs.sh`     Predict optimal similarity cut-offs
  --------------------------------------------------------------------------------

------------------------------------------------------------------------

### R modules

  Script                    Description
  ------------------------- ----------------------------------------------------
  `R/utils.R`               Shared utility functions (e.g., `is_identified()`)
  `R/reformat.R`            FASTA header parsing and taxonomy extraction
  `R/check_annotations.R`   Infraspecific annotation standardisation
  `R/subset.R`              Balanced taxon subset generation
  `R/compute_sim.R`         Pairwise similarity computation using vsearch
  `R/predict.R`             Cut-off prediction (parallel or sequential)

------------------------------------------------------------------------

## Directory structure

    data/
    ├── full_ITS/    # Full ITS sequences, taxonomy, ID files, predictions
    ├── ITS1/        # ITS1 sequences, taxonomy, ID files, predictions
    └── ITS2/        # ITS2 sequences, taxonomy, ID files, predictions

------------------------------------------------------------------------

## Environment setup

A conda environment specification is provided in `environment.yml`.

[mamba](https://mamba.readthedocs.io/) is strongly recommended over
`conda` for faster dependency resolution.

### Local installation

``` bash
mamba env create -f environment.yml
conda activate dyna_clust_predict
```

### HPC (SLURM)

Run from the project root:

``` bash
sbatch scripts/create_dyna_clust_env.sh
conda activate dyna_clust_predict
```

------------------------------------------------------------------------

## Quick start

From the project root:

``` bash
conda activate dyna_clust_predict

sbatch scripts/01_reformat_ITS.sh
sbatch scripts/02_check_annotations.sh
sbatch scripts/03_extract_subregions.sh
sbatch scripts/04_prepare_subsets.sh
sbatch scripts/05_compute_sim.sh   # Optional
sbatch scripts/06_predict_cut-offs.sh
```

> **Note:** Step 05 is optional. Similarity can be computed on-the-fly
> in Step 06, which is preferred for large datasets and parallel
> execution.

------------------------------------------------------------------------

# Key parameters

## Sequence selection

Used in: `subset.R` (via `04_prepare_subsets.sh`)

These parameters control taxonomic balance and sampling constraints
prior to similarity prediction.

    --min_subgroups   INT    Minimum unique child taxa per parent (default: 10)
    --min_sequences   INT    Minimum sequences per parent after proportion cap (default: 30)
    --max_sequences   INT    Maximum sequences per parent; excess downsampled proportionally (default: 25000)
    --max_proportion  FLOAT  Maximum fraction a child taxon may represent (default: 0.5)

Example:

``` bash
Rscript R/subset.R \
  --min_subgroups 15 \
  --min_sequences 50 \
  --max_sequences 20000 \
  --max_proportion 0.4
```

------------------------------------------------------------------------

## Similarity prediction

Used in: `predict.R` (via `06_predict_cut-offs.sh`)

The same parameters are applied during cut-off optimisation to ensure
taxonomic balance during F-measure calculation.

    --min_subgroups
    --min_sequences
    --max_sequences
    --max_proportion

Defaults are identical to those used during subset preparation.

------------------------------------------------------------------------

## Citation

Vu, D., Nilsson, R. H., & Verkley, G. J. M. (2022).\
*dnabarcoder: an open-source software package for analyzing and
predicting DNA sequence similarity cut-offs for fungal sequence
identification.*\
Molecular Ecology Resources.\
https://doi.org/10.1111/1755-0998.13651
