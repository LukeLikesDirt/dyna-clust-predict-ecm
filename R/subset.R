#!/usr/bin/env Rscript
# subset.R — Prepare sequence subsets and ID files for prediction
#
# Only unique sequences identified to species rank are used as the 
# base pool for all ID-list generation steps
#
# STEP 1 — Nested prediction ID lists:
#   For each (target x parent) combination, apply sequential filters:
#     1. min_subgroups  — min resolved child taxa per parent chunk
#     2. max_proportion — cap dominant child taxon via random subsampling
#     3. min_sequences  — min sequences after capping
#     4. max_sequences  — random downsample if too large
#
# STEP 2 — Global prediction ID lists:
#   For each target rank, cap the dominant clade via random subsampling,
#   then random downsample to max_sequences.
#
# Usage:
#   Rscript subset.R --classification_in taxonomy.tsv \
#                    --output_dir data/full_ITS
#
# Note: This script must be run from the project root directory.

required_packages <- c("optparse", "readr", "dplyr", "data.table")
missing_packages  <- required_packages[
  !sapply(required_packages, requireNamespace, quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Missing packages: ", paste(missing_packages, collapse = ", "),
       "\nInstall with: install.packages(c('", paste(missing_packages, collapse = "', '"), "'))")
}

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(data.table)
})

source("R/utils.R")

# ── Arguments ─────────────────────────────────────────────────────────────────

option_list <- list(
  make_option("--classification_in",
              type = "character", default = "data/full_ITS/eukaryome_ITS.classification",
              metavar = "FILE",
              help = "Input tab-delimited classification file [default: %default]"),
  make_option("--output_dir",
              type = "character", default = "data/full_ITS",
              metavar = "DIR",
              help = "Directory for output ID files [default: %default]"),
  make_option("--min_subgroups",
              type = "integer", default = 10L, metavar = "INT",
              help = "Min unique child taxa per parent chunk [default: %default]"),
  make_option("--min_sequences",
              type = "integer", default = 30L, metavar = "INT",
              help = "Min sequences per chunk after proportion cap [default: %default]"),
  make_option("--max_sequences",
              type = "integer", default = 20000L, metavar = "INT",
              help = "Max sequences per chunk; excess is randomly downsampled [default: %default]"),
  make_option("--max_proportion",
              type = "double", default = 1, metavar = "NUM",
              help = "Max fraction of chunk that dominant child taxon may represent [default: %default]")
)

opt <- parse_args(
  OptionParser(
    option_list = option_list,
    usage       = "%prog [options]",
    description = paste(
      "Generate prediction ID files for dyna-clust-predict.",
      "Produces nested and global ID files for each rank."
    )
  )
)

classification_in <- opt$classification_in
output_dir        <- opt$output_dir
min_subgroups     <- opt$min_subgroups
min_sequences     <- opt$min_sequences
max_sequences     <- opt$max_sequences
max_proportion    <- opt$max_proportion

if (!file.exists(classification_in)) stop("Classification file not found: ", classification_in)
if (!is.numeric(min_subgroups) || min_subgroups <= 0)
  stop("min_subgroups must be a positive integer")
if (!is.numeric(min_sequences) || min_sequences <= 0)
  stop("min_sequences must be a positive integer")
if (!is.numeric(max_sequences) || max_sequences <= 0)
  stop("max_sequences must be a positive integer")
if (is.na(max_proportion) || max_proportion <= 0 || max_proportion >= 1)
  stop("max_proportion must be strictly between 0 and 1")

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

cat("=== PARAMETERS ===\n")
cat(sprintf("  classification_in : %s\n", classification_in))
cat(sprintf("  output_dir        : %s\n", output_dir))
cat(sprintf("  min_subgroups     : %d\n", min_subgroups))
cat(sprintf("  min_sequences     : %d\n", min_sequences))
cat(sprintf("  max_sequences     : %d\n", max_sequences))
cat(sprintf("  max_proportion    : %.2f\n", max_proportion))
cat("\n")

# ── Rank metadata ─────────────────────────────────────────────────────────────

rank_hierarchy <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")

rank_abbr <- c(
  kingdom = "kng", phylum = "phy", class = "cls", order = "ord",
  family  = "fam", genus  = "gen", species = "spe"
)

parent_ranks_map <- list(
  species = c("genus", "family", "order", "class", "phylum", "kingdom"),
  genus   = c("family", "order", "class", "phylum", "kingdom"),
  family  = c("order", "class", "phylum", "kingdom"),
  order   = c("class", "phylum", "kingdom"),
  class   = c("phylum", "kingdom"),
  phylum  = c("kingdom"),
  kingdom = character(0)
)

# ── Functions ─────────────────────────────────────────────────────────────────

# nested_prediction_filter: apply sequential filters per parent chunk.
nested_prediction_filter <- function(df, target_rank, parent_rank,
                                     max_proportion = 1, min_subgroups = 10,
                                     min_sequences = 30, max_sequences = 20000) {
  df <- df %>%
    filter(is_identified(!!sym(parent_rank)), is_identified(!!sym(target_rank)))
  if (nrow(df) == 0) return(NULL)

  df_split      <- split(df, df[[parent_rank]])
  chunk_results <- list()

  for (parent_taxon in names(df_split)) {
    chunk <- df_split[[parent_taxon]]

    # Filter 1: min unique child taxa
    if (length(unique(chunk[[target_rank]])) < min_subgroups) next

    # Filter 2: cap dominant child taxon via random subsampling
    child_counts   <- chunk %>%
      group_by(!!sym(target_rank)) %>%
      summarise(n = n(), .groups = "drop") %>%
      arrange(desc(n))

    dominant_child <- child_counts[[1, target_rank]]
    dominant_prop  <- child_counts[[1, "n"]] / nrow(chunk)

    if (dominant_prop > max_proportion) {
      non_dom <- chunk %>% filter(!!sym(target_rank) != dominant_child)
      dom     <- chunk %>% filter(!!sym(target_rank) == dominant_child)
      max_dom <- floor(nrow(non_dom) * max_proportion / (1 - max_proportion))
      dom     <- slice_sample(dom, n = max_dom)
      chunk   <- bind_rows(non_dom, dom)
    }

    # Filter 3: min sequences after capping
    if (nrow(chunk) < min_sequences) next

    # Filter 4: random downsample if too large
    if (nrow(chunk) > max_sequences) chunk <- slice_sample(chunk, n = max_sequences)

    chunk_results[[parent_taxon]] <- chunk
  }

  if (length(chunk_results) == 0) return(NULL)
  bind_rows(chunk_results)
}

# global_prediction_filter: cap the dominant target-rank clade via random
# subsampling, then randomly downsample to max_sequences.
global_prediction_filter <- function(df, target_rank, max_proportion = 1,
                                     max_sequences = 20000, min_sequences = 30) {
  df <- df %>% filter(is_identified(!!sym(target_rank)))
  if (nrow(df) < min_sequences) return(NULL)

  clade_counts   <- df %>%
    group_by(!!sym(target_rank)) %>%
    summarise(n = n(), .groups = "drop") %>%
    arrange(desc(n))

  dominant_clade <- clade_counts[[1, target_rank]]
  dominant_prop  <- clade_counts[[1, "n"]] / nrow(df)

  if (dominant_prop > max_proportion) {
    non_dom <- df %>% filter(!!sym(target_rank) != dominant_clade)
    dom     <- df %>% filter(!!sym(target_rank) == dominant_clade)
    max_dom <- floor(nrow(non_dom) * max_proportion / (1 - max_proportion))
    dom     <- slice_sample(dom, n = max_dom)
    df      <- bind_rows(non_dom, dom)
  }

  if (nrow(df) < min_sequences) return(NULL)
  if (nrow(df) > max_sequences) df <- slice_sample(df, n = max_sequences)
  if (nrow(df) < min_sequences) return(NULL)
  df
}

# ── Read inputs ───────────────────────────────────────────────────────────────

cat("Reading classification file...\n")
classification_df <- fread(classification_in) %>%
  select(id, kingdom, phylum, class, order, family, genus, species)

# ── Pre-filter: species-identified  ───────────────────────────────────────────

cat("\nPre-filtering to species-identified sequences (>= 3 per species)...\n")
n_before <- nrow(classification_df)

classification_df <- classification_df %>%
  filter(is_identified(species))

cat(sprintf("  Retained %d / %d sequences (%d species)\n\n",
            nrow(classification_df), n_before,
            length(unique(classification_df$species))))

# ── STEP 1: Nested prediction ID lists ───────────────────────────────────────

cat("STEP 1: Preparing nested prediction ID lists...\n")

for (target_rank in names(parent_ranks_map)) {
  valid_parents <- parent_ranks_map[[target_rank]]
  if (length(valid_parents) == 0) {
    cat(sprintf("  %-10s  no valid parent ranks — skipping\n", target_rank)); next
  }
  for (parent_rank in valid_parents) {
    out_path <- file.path(output_dir,
                          sprintf("%s_pred_id_%s.txt", target_rank, rank_abbr[[parent_rank]]))
    result <- nested_prediction_filter(
      df             = classification_df, target_rank = target_rank, parent_rank = parent_rank,
      max_proportion = max_proportion, min_subgroups = min_subgroups,
      min_sequences  = min_sequences,  max_sequences  = max_sequences
    )
    if (is.null(result)) {
      cat(sprintf("  %-10s within %-10s ->  no groups passed filters\n",
                  target_rank, parent_rank)); next
    }
    writeLines(result$id, out_path)
    cat(sprintf("  %-10s within %-10s ->  %6d IDs  |  %d unique %s  ->  %s\n",
                target_rank, parent_rank, nrow(result),
                length(unique(result[[target_rank]])), target_rank, out_path))
  }
}

# ── STEP 2: Global prediction ID lists ───────────────────────────────────────

cat("\nSTEP 2: Preparing global prediction ID lists...\n")

for (target_rank in rank_hierarchy) {
  out_path <- file.path(output_dir, sprintf("%s_pred_id_global.txt", target_rank))
  result   <- global_prediction_filter(
    df = classification_df, target_rank = target_rank,
    max_proportion = max_proportion, max_sequences = max_sequences, min_sequences = min_sequences
  )
  if (is.null(result)) {
    cat(sprintf("  %-10s ->  no sequences passed filters\n", target_rank)); next
  }
  writeLines(result$id, out_path)
  dom_pct <- max(table(result[[target_rank]])) / nrow(result) * 100
  cat(sprintf(
    "  %-10s ->  %6d IDs  |  %d unique %-10s  |  dominant: %.1f%%  ->  %s\n",
    target_rank, nrow(result), length(unique(result[[target_rank]])),
    target_rank, dom_pct, out_path
  ))
}

cat("\nDone. All ID files saved to:", output_dir, "\n")
