# dyna-clust-predict-ecm

`dyna-clust-predict-ecm` is a **special-use adaptation** of

`dyna-clust-predict` for predicting similarity cut-offs in
**ectomycorrhizal (ECM) fungi** using the
**full-ITS, ITS or ITS2** fragments.

It predicts optimal sequence similarity cut-offs for classification and
clustering using **vsearch global alignment** and
**F-measure optimisation**, following Vu *et al.* (2022).

## Overview

`dyna-clust-predict-ecm` predicts **species-level** similarity cut-offs
for ectomycorrhizal (ECM) fungi using a curated ITS reference database
(EUKARYOME). ECM fungi are identified based on genera listed in
FungalTraits (Põlme *et al.* 2020), and lineage-level similarity cut-offs are
leveraged from Tedersoo *et al.* (2022).

This is an early release and is currently only configured for full-ITS
sequences. The pipeline can be configured to support ITS1 and ITS2
subregions as well, and these regions will be added in a future release.

The pipeline:

- Extracts ITS1 and ITS2 subregions using ITSx
- Computes global pairwise similarity using vsearch
- Optimises similarity cut-offs using F-measure
- Supports parallel processing for efficient large-scale analyses

The workflow is conceptually adapted from
[dnabarcoder](https://github.com/vuthuyduong/dnabarcoder), but:

- Uses **vsearch** for global alignment
- Supports scalable parallel computation
- Integrates similarity prediction within an R-based workflow

## Pipeline steps

  `scripts/01_reformat_ITS.sh`        	Download EUKARYOME, reformat headers,
                                       	extract taxonomy

  `scripts/02_check_annotations.sh`    	Standardise infraspecific annotations

  `scripts/03_extract_subregions.sh`   	Extract ITS1 and ITS2 using ITSx

  `scripts/04_subset_ecm.R`            	Filter ECM taxa and generate balanced
                                       	prediction subsets

  `scripts/05_compute_sim.sh`          	*(Optional)* Pre-compute similarity
                                       	matrices

  `scripts/06_predict_ecm_cutoffs.sh`  Predict optimal similarity cut-offs

All scripts must be run from the **project root directory**.

### R modules

  `R/utils.R`             Shared utility functions (e.g., `is_identified()`)

  `R/reformat.R`          FASTA header parsing and taxonomy extraction

  `R/check_annotations.R`	Infraspecific annotation standardisation

  `R/subset.R`            Balanced taxon subset generation

  `R/compute_sim.R`       Pairwise similarity computation using vsearch

  `R/predict.R`           Cut-off prediction (parallel or sequential)

## Citations

Põlme, S., Abarenkov, K., Henrik Nilsson, R., Lindahl, B. D., Clemmensen, K. E.,
Kauserud, H., ... & Tedersoo, L. (2020). FungalTraits: a user-friendly traits
database of fungi and fungus-like stramenopiles. *Fungal Diversity*, *105*(1),
1–16. [https://doi.org/10.1007/s13225-020-00466-2](https://doi.org/10.1007/s13225-020-00466-2)

Tedersoo, L., Bahram, M., Zinger, L., Nilsson, R. H., Kennedy, P. G., Yang, T.,
Anslan, S. & Mikryukov, V. (2022). Best practices in metabarcoding of fungi: from
experimental design to results. *Molecular Ecology*, *31*(10), 2769–2795.

Tedersoo, L., Hosseyni Moghaddam, M. S., Mikryukov, V., Hakimzadeh, A., Bahram,
M., Nilsson, R. H., ... & Anslan, S. (2024). EUKARYOME: the rRNA gene reference
database for identification of all eukaryotes. Database, 2024, baae043.
https://doi.org/10.1093/database/baae043

Vu, D., Nilsson, R. H., & Verkley, G. J. (2022). Dnabarcoder: An open‐source
software package for analysing and predicting DNA sequence similarity cutoffs for
fungal sequence identification. *Molecular Ecology Resources*, *22*(7), 2793–2809.
[https://doi.org/10.1111/1755-0998.13651](https://doi.org/10.1111/1755-0998.13651)
