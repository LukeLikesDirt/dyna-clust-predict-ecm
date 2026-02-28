#!/usr/bin/env Rscript
# subset.R — Prepare sequence subsets and ID files for prediction
#
# Produces three sets of plain-text ID files (one sequence ID per line) so
# the global FASTA and classification table can be filtered at prediction time.
#
# STEP 1 — Unique sequences:
#   Remove duplicate sequences at each rank, keeping one representative per
#   taxonomic group. Outputs one ID file per rank.
#
# STEP 2 — Nested prediction ID lists:
#   For each (target x parent) combination, apply four sequential filters:
#     1. min_subgroups  — min resolved child taxa per parent chunk
#     2. max_proportion — cap dominant child taxon to <= max_proportion
#     3. min_sequences  — min sequences after capping
#     4. max_sequences  — balanced round-robin downsample if too large
#
# STEP 3 — Global prediction ID lists:
#   Top-down recursive budget allocation across the full hierarchy, ensuring
#   rare clades are represented proportionally at every taxonomic level.
#
# Usage:
#   Rscript subset.R --fasta_in sequences.fasta \
#                    --classification_in taxonomy.tsv \
#                    --output_dir data/full_ITS
#
# Note: This script must be run from the project root directory.

required_packages <- c("optparse", "readr", "Biostrings", "dplyr", "data.table")
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
  library(Biostrings)
  library(dplyr)
  library(data.table)
})

source("R/utils.R")

# ── Arguments ─────────────────────────────────────────────────────────────────

option_list <- list(
  make_option("--fasta_in",
              type = "character", default = "data/full_ITS/eukaryome_ITS.fasta",
              metavar = "FILE",
              help = "Input FASTA file [default: %default]"),
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
              type = "integer", default = 25000L, metavar = "INT",
              help = "Max sequences per chunk; excess is balanced-downsampled [default: %default]"),
  make_option("--max_proportion",
              type = "double", default = 0.5, metavar = "NUM",
              help = "Max fraction of chunk that dominant child taxon may represent [default: %default]")
)

opt <- parse_args(
  OptionParser(
    option_list = option_list,
    usage       = "%prog [options]",
    description = paste(
      "Generate prediction ID files for dyna-clust-predict.",
      "Produces unique-sequence, nested, and global ID files for each rank."
    )
  )
)

fasta_in          <- opt$fasta_in
classification_in <- opt$classification_in
output_dir        <- opt$output_dir
min_subgroups     <- opt$min_subgroups
min_sequences     <- opt$min_sequences
max_sequences     <- opt$max_sequences
max_proportion    <- opt$max_proportion

if (!file.exists(fasta_in))          stop("FASTA file not found: ", fasta_in)
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
cat(sprintf("  fasta_in          : %s\n", fasta_in))
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

# filter_unique_sequences: removes duplicate sequences at a given rank,
# keeping the representative from the most abundant taxon per duplicate group.
filter_unique_sequences <- function(classification_df, fasta_seqs, rank = "species") {

  valid_ranks <- rank_hierarchy
  if (!rank %in% valid_ranks) stop("Invalid rank: ", rank)

  rank_hier <- valid_ranks[seq_len(which(valid_ranks == rank))]

  rank_df   <- classification_df %>% filter(is_identified(.data[[rank]]))
  group_var <- if (rank == "kingdom") "kingdom" else rank_hier[length(rank_hier) - 1]
  rank_df   <- rank_df %>% group_by(across(all_of(group_var)))

  rank_seqs <- fasta_seqs[names(fasta_seqs) %in% rank_df$id]

  duplicated_seqs <- tibble(id = names(rank_seqs), sequence = as.character(rank_seqs)) %>%
    group_by(sequence) %>%
    filter(n() > 1) %>%
    mutate(group_id = paste0("group_", cur_group_id())) %>%
    ungroup() %>%
    select(-sequence)

  if (nrow(duplicated_seqs) == 0) {
    message("No duplicated sequences at ", rank, " level")
    result <- list(rank_df, rank_seqs)
    names(result) <- c(paste0(rank, "_df"), paste0(rank, "_seqs"))
    return(result)
  }

  unique_seqs <- duplicated_seqs %>%
    left_join(classification_df %>% select(id, all_of(rank_hier)), by = "id") %>%
    group_by(across(c(group_id, all_of(rank_hier)))) %>%
    mutate(n_seqs = n()) %>%
    filter(n_seqs == max(n_seqs)) %>%
    arrange(desc(n_seqs)) %>%
    ungroup() %>%
    group_by(group_id) %>%
    slice(1)

  seqs_to_remove <- setdiff(duplicated_seqs$id, unique_seqs$id)
  rank_df        <- rank_df %>% filter(!id %in% seqs_to_remove)
  rank_seqs      <- rank_seqs[names(rank_seqs) %in% rank_df$id]

  result <- list(rank_df, rank_seqs)
  names(result) <- c(paste0(rank, "_df"), paste0(rank, "_seqs"))
  result
}

# balanced_downsample: round-robin downsample to n rows, balanced across taxa.
balanced_downsample <- function(df, target_rank, n) {
  df %>%
    ungroup() %>%
    group_by(!!sym(target_rank)) %>%
    mutate(.draw_order = row_number()) %>%
    ungroup() %>%
    arrange(.draw_order, !!sym(target_rank)) %>%
    slice_head(n = n) %>%
    select(-.draw_order)
}

# nested_prediction_filter: apply four sequential filters per parent chunk.
nested_prediction_filter <- function(df, target_rank, parent_rank,
                                     max_proportion = 0.5, min_subgroups = 10,
                                     min_sequences = 30, max_sequences = 25000) {
  df <- df %>%
    filter(is_identified(!!sym(parent_rank)), is_identified(!!sym(target_rank)))
  if (nrow(df) == 0) return(NULL)

  df_split      <- split(df, df[[parent_rank]])
  chunk_results <- list()

  for (parent_taxon in names(df_split)) {
    chunk <- df_split[[parent_taxon]]

    # Filter 1: min unique child taxa
    if (length(unique(chunk[[target_rank]])) < min_subgroups) next

    # Filter 2: cap dominant child taxon to max_proportion
    rank_counts <- chunk %>%
      group_by(!!sym(target_rank)) %>%
      summarise(n = n(), .groups = "drop") %>%
      arrange(desc(n))

    if (nrow(rank_counts) == 1) {
      selected_ids <- chunk$id
    } else {
      largest_taxon    <- rank_counts[[1, target_rank]]
      smaller_count    <- sum(rank_counts$n[-1])
      max_from_largest <- floor((max_proportion / (1 - max_proportion)) * smaller_count)
      non_dom_ids <- chunk %>% filter(!!sym(target_rank) != largest_taxon) %>% pull(id)
      dom_ids     <- chunk %>% filter(!!sym(target_rank) == largest_taxon) %>% pull(id)
      selected_ids <- if (length(dom_ids) <= max_from_largest) {
        c(non_dom_ids, dom_ids)
      } else {
        c(non_dom_ids, sample(dom_ids, max_from_largest))
      }
    }
    chunk <- chunk %>% filter(id %in% selected_ids)

    # Filter 3: min sequences after capping
    if (nrow(chunk) < min_sequences) next

    # Filter 4: balanced downsample if too large
    if (nrow(chunk) > max_sequences) chunk <- balanced_downsample(chunk, target_rank, max_sequences)

    chunk_results[[parent_taxon]] <- chunk
  }

  if (length(chunk_results) == 0) return(NULL)
  bind_rows(chunk_results)
}

# allocate_budget: distribute budget across groups, redistributing from
# underfull groups to groups that can absorb more.
allocate_budget <- function(available, budget) {
  n <- length(available)
  if (n == 0 || budget == 0) return(setNames(integer(n), names(available)))
  if (sum(available) <= budget) return(available)

  allocation <- setNames(integer(n), names(available))
  remaining  <- budget
  competing  <- rep(TRUE, n)

  while (remaining > 0 && any(competing)) {
    n_competing <- sum(competing)
    per_group   <- remaining / n_competing
    underfull   <- competing & (available < per_group)

    if (!any(underfull)) {
      base_share              <- floor(remaining / n_competing)
      allocation[competing]  <- base_share
      leftover               <- remaining - base_share * n_competing
      if (leftover > 0) {
        top_up_idx              <- which(competing)[seq_len(leftover)]
        allocation[top_up_idx] <- allocation[top_up_idx] + 1L
      }
      break
    }
    for (i in which(underfull)) {
      allocation[i] <- available[i]
      remaining     <- remaining - available[i]
      competing[i]  <- FALSE
    }
  }
  allocation
}

# hierarchical_sample: recursive budget descent through taxonomy.
hierarchical_sample <- function(df, grouping_ranks, target_rank, budget) {
  if (nrow(df) == 0 || budget == 0) return(NULL)
  if (length(grouping_ranks) == 0) {
    return(balanced_downsample(df, target_rank, min(budget, nrow(df))))
  }

  current_rank    <- grouping_ranks[1]
  remaining_ranks <- grouping_ranks[-1]
  df              <- df %>% filter(is_identified(!!sym(current_rank)))
  if (nrow(df) == 0) return(NULL)

  groups      <- split(df, df[[current_rank]])
  available   <- sapply(groups, nrow)
  allocations <- allocate_budget(available, budget)

  results <- list()
  for (grp in names(groups)) {
    if (allocations[[grp]] == 0) next
    res <- hierarchical_sample(groups[[grp]], remaining_ranks, target_rank, allocations[[grp]])
    if (!is.null(res) && nrow(res) > 0) results[[grp]] <- res
  }
  if (length(results) == 0) return(NULL)
  bind_rows(results)
}

# global_prediction_filter: top-down recursive budget allocation.
global_prediction_filter <- function(df, target_rank, hierarchy = rank_hierarchy,
                                     max_sequences = 25000, min_sequences = 30) {
  target_idx     <- which(hierarchy == target_rank)
  grouping_ranks <- hierarchy[seq_len(target_idx - 1)]

  df <- df %>% filter(is_identified(!!sym(target_rank)))
  if (nrow(df) < min_sequences) return(NULL)

  if (length(grouping_ranks) == 0) {
    result <- balanced_downsample(df, target_rank, min(max_sequences, nrow(df)))
    if (nrow(result) < min_sequences) return(NULL)
    return(result)
  }

  result <- hierarchical_sample(df, grouping_ranks, target_rank, max_sequences)
  if (is.null(result) || nrow(result) < min_sequences) return(NULL)
  result
}

# ── Read inputs ───────────────────────────────────────────────────────────────

cat("Reading classification file...\n")
classification_df <- fread(classification_in) %>%
  select(id, kingdom, phylum, class, order, family, genus, species)

cat("Reading FASTA file...\n")
fasta_seqs <- readDNAStringSet(fasta_in)
cat(sprintf("  Sequences loaded: %d\n\n", length(fasta_seqs)))

# ── STEP 1: Unique sequences ──────────────────────────────────────────────────

cat("STEP 1: Selecting unique sequences...\n")

rank_dfs <- list()
for (rank in rank_hierarchy) {
  result           <- filter_unique_sequences(classification_df, fasta_seqs, rank = rank)
  rank_dfs[[rank]] <- result[[paste0(rank, "_df")]]
  out_path         <- file.path(output_dir, sprintf("%s_unique_id.txt", rank))
  writeLines(rank_dfs[[rank]]$id, out_path)
  cat(sprintf("  %-10s  %6d unique sequences  ->  %s\n",
              rank, nrow(rank_dfs[[rank]]), out_path))
}

# ── STEP 2: Nested prediction ID lists ───────────────────────────────────────

cat("\nSTEP 2: Preparing nested prediction ID lists...\n")

for (target_rank in names(parent_ranks_map)) {
  valid_parents <- parent_ranks_map[[target_rank]]
  if (length(valid_parents) == 0) {
    cat(sprintf("  %-10s  no valid parent ranks — skipping\n", target_rank)); next
  }
  df_target <- rank_dfs[[target_rank]]
  for (parent_rank in valid_parents) {
    out_path <- file.path(output_dir,
                          sprintf("%s_pred_id_%s.txt", target_rank, rank_abbr[[parent_rank]]))
    result <- nested_prediction_filter(
      df             = df_target, target_rank = target_rank, parent_rank = parent_rank,
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

# ── STEP 3: Global prediction ID lists ───────────────────────────────────────

cat("\nSTEP 3: Preparing global prediction ID lists (hierarchical budget allocation)...\n")

for (target_rank in rank_hierarchy) {
  out_path <- file.path(output_dir, sprintf("%s_pred_id_global.txt", target_rank))
  result   <- global_prediction_filter(
    df = rank_dfs[[target_rank]], target_rank = target_rank,
    hierarchy = rank_hierarchy, max_sequences = max_sequences, min_sequences = min_sequences
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
