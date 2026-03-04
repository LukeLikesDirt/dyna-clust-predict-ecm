#!/usr/bin/env Rscript
# predict.R — Predict optimal similarity cutoffs for DNA barcoding taxonomy.
#
# This scritp is adapted from  the predict function in dnabarcoder: 
# https://github.com/vuthuyduong/dnabarcoder
#
# Uses vsearch --allpairs_global (global alignment) for similarity computation.
# A pre-computed .sim file may be supplied via --sim to skip that step.
#
# Supports both sequential and parallel (furrr/future) dataset processing,
# controlled by --run_parallel (default: yes).
#
#
# Usage:
#   Rscript R/predict.R \
#     --input sequences.fasta \
#     --classification taxonomy.tsv \
#     --rank species \
#     --start_threshold 0.9 \
#     --end_threshold 1.0 \
#     --run_parallel yes \
#     --n_cpus 8
#
# Algorithm:
#   1. For each taxonomic rank (e.g. species):
#      a. Partition sequences into datasets (one global, or one per higher-rank
#         taxon when --higher_rank is set).
#      b. Sweep similarity thresholds from --start_threshold to --end_threshold.
#      c. At each threshold build a neighbour graph (pairs with sim >= threshold).
#      d. Find connected components (predicted clusters) via iterative BFS.
#      e. Compute F-measure between predicted clusters and true taxonomic classes.
#      f. Select the threshold that maximises the F-measure.
#   2. Save full prediction traces to JSON and optimal cutoffs to JSON + TSV.
#
# Outputs (in --out directory):
#   <prefix>.predicted         — full F-measure traces per threshold (JSON)
#   <prefix>.cutoffs.json      — optimal cutoffs only (JSON)
#   <prefix>.cutoffs.json.txt  — tab-delimited cutoff summary
#
# Required external tool: vsearch (https://github.com/torognes/vsearch)

# ── Packages ──────────────────────────────────────────────────────────────────
# All packages are always loaded regardless of --run_parallel.

suppressPackageStartupMessages({
  library(optparse)
  library(jsonlite)
  library(data.table)
  library(furrr)
  library(future)
})

# Allow large globals to be exported to workers (sim_dt can be several hundred MB).
options(future.globals.maxSize = 2000 * 1024^2)

# ── Shared utilities ──────────────────────────────────────────────────────────
source("R/utils.R")

# ── Arguments ─────────────────────────────────────────────────────────────────

option_list <- list(
  make_option(c("-i", "--input"),
              type = "character", metavar = "FILE",
              help = "Input FASTA file [required]"),
  make_option(c("-c", "--classification"),
              type = "character", default = "", metavar = "FILE",
              help = "Tab-delimited classification file [required]"),
  make_option(c("-o", "--out"),
              type = "character", default = "dnabarcoder", metavar = "DIR",
              help = "Output directory [default: %default]"),
  make_option("--prefix",
              type = "character", default = "", metavar = "STR",
              help = "Prefix for output filenames [default: FASTA basename]"),
  make_option("--rank",
              type = "character", default = "", metavar = "STR",
              help = "Target rank(s) to predict, comma-separated (e.g. species,genus) [required]"),
  make_option("--higher_rank",
              type = "character", default = "", metavar = "STR",
              help = paste("Higher rank(s) for local prediction, comma-separated (e.g. genus,family).",
                           "Omit for a single global prediction across all sequences.",
                           "[default: \"\"]")),
  make_option("--start_threshold",
              type = "double", default = 0, metavar = "NUM",
              help = "Starting similarity threshold [default: %default]"),
  make_option("--end_threshold",
              type = "double", default = 1, metavar = "NUM",
              help = "Ending similarity threshold [default: %default]"),
  make_option("--step",
              type = "double", default = 0.001, metavar = "NUM",
              help = "Threshold step size [default: %default]"),
  make_option("--sim",
              type = "character", default = "", metavar = "FILE",
              help = "Pre-computed .sim file. Computed with vsearch if not supplied. [default: \"\"]"),
  make_option("--min_group_no",
              type = "integer", default = 10L, metavar = "INT",
              help = "Min taxonomic groups required to report a cutoff [default: %default]"),
  make_option("--min_seq_no",
              type = "integer", default = 30L, metavar = "INT",
              help = "Min sequences required to report a cutoff [default: %default]"),
  make_option("--max_seq_no",
              type = "integer", default = 25000L, metavar = "INT",
              help = "Max sequences per dataset; excess is randomly sampled [default: %default]"),
  make_option("--max_proportion",
              type = "double", default = 1.0, metavar = "NUM",
              help = paste("Skip datasets where the dominant group exceeds this proportion",
                           "[default: %default]")),
  make_option("--min_cutoff",
              type = "double", default = 0, metavar = "NUM",
              help = "Min cutoff value to include in output [default: %default]"),
  make_option("--id_col",
              type = "character", default = "id", metavar = "STR",
              help = "ID column name in the classification file [default: %default]"),
  make_option("--n_cpus",
              type = "integer",
              default = max(1L, future::availableCores() - 1L),
              metavar = "INT",
              help = paste("Workers (parallel) or vsearch threads (sequential).",
                           "[default: availableCores() - 1]")),
  make_option("--run_parallel",
              type = "character", default = "yes", metavar = "yes/no",
              help = "Run datasets in parallel using furrr/future [default: %default]"),
  make_option("--redo",
              type = "character", default = "no", metavar = "yes/no",
              help = "Recompute F-measures even if cached [default: %default]"),
  make_option("--tmp_dir",
              type = "character", default = "./tmp", metavar = "DIR",
              help = "Directory for temporary vsearch output [default: %default]")
)

opt <- parse_args(
  OptionParser(
    option_list = option_list,
    usage       = "%prog [options] -i sequences.fasta -c taxonomy.tsv --rank species",
    description = paste(
      "Predict optimal similarity cutoffs for DNA barcoding taxonomy.",
      "Uses vsearch --allpairs_global (global alignment) for similarity computation.",
      "Sweeps thresholds, clusters sequences, computes F-measure vs known taxonomy,",
      "and selects the threshold maximising F-measure.",
      "Parallel or sequential dataset processing is controlled by --run_parallel."
    )
  )
)

# ── Boolean flags from string options ────────────────────────────────────────

run_parallel <- identical(tolower(opt$run_parallel), "yes")
redo         <- identical(tolower(opt$redo),         "yes")

# ── Validation ────────────────────────────────────────────────────────────────

if (is.null(opt$input))           stop("--input is required.")
if (!file.exists(opt$input))      stop("FASTA file not found: ", opt$input)
if (opt$rank == "")               stop("--rank is required (e.g. --rank species).")
if (nchar(opt$classification) == 0) stop("--classification is required.")

# ── Assign parameters to short names ─────────────────────────────────────────

fasta_file      <- opt$input
class_file      <- opt$classification
output_dir      <- opt$out
prefix          <- opt$prefix
rank_str        <- opt$rank
higher_rank_str <- opt$higher_rank
start_t         <- opt$start_threshold
end_t           <- opt$end_threshold
step_t          <- opt$step
sim_file_arg    <- opt$sim
min_group_no    <- opt$min_group_no
min_seq_no      <- opt$min_seq_no
max_seq_no      <- opt$max_seq_no
max_prop_limit  <- opt$max_proportion
min_cutoff      <- opt$min_cutoff
id_col          <- opt$id_col
n_cpus          <- opt$n_cpus
tmp_dir         <- opt$tmp_dir

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

if (prefix == "" || is.null(prefix)) {
  prefix <- sub("\\.[^./]+$", "", basename(fasta_file))
}
prediction_file <- file.path(output_dir, paste0(prefix, ".predicted"))

# ── Parallel backend setup ─────────────────────────────────────────────────────
# multicore (fork) on Unix/macOS — shares memory, no serialisation overhead.
# multisession on Windows — serialises globals, works everywhere.
# Sequential mode leaves vsearch threads available to the single worker.

if (run_parallel) {
  if (.Platform$OS.type == "unix") {
    plan(multicore,    workers = n_cpus)
  } else {
    plan(multisession, workers = n_cpus)
  }
  cat(sprintf("[predict] Parallel backend: %s  workers: %d\n",
              if (.Platform$OS.type == "unix") "multicore" else "multisession", n_cpus))
} else {
  plan(sequential)
  cat(sprintf("[predict] Sequential mode  vsearch threads: %d\n", n_cpus))
}

# ── NULL-coalescing operator ──────────────────────────────────────────────────

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

# ══════════════════════════════════════════════════════════════════════════════
# Functions
# ══════════════════════════════════════════════════════════════════════════════

# ── Load FASTA — extract sequence IDs and header descriptions ─────────────────

load_fasta_ids <- function(fasta_file) {
  lines <- readLines(fasta_file)
  h_idx <- which(startsWith(lines, ">"))
  if (length(h_idx) == 0) stop("No sequences found in FASTA file: ", fasta_file)
  ids   <- sub("^>([^ ]+).*", "\\1", lines[h_idx])
  descs <- sub("^>", "",            lines[h_idx])
  names(descs) <- ids
  list(ids = ids, descriptions = descs)
}

# ── Load similarity matrix from a .sim file ───────────────────────────────────
# Format: id1 id2 score (space-delimited, no header).
# The file stores only one direction per pair; this function symmetrises the
# matrix and adds self-similarities (score = 1.0) so downstream code can
# filter by seq_ids without special-casing the diagonal.

load_sim <- function(sim_file) {
  cat("[predict] Loading similarity matrix from:", sim_file, "\n")
  sim_dt <- fread(sim_file, header = FALSE, col.names = c("i", "j", "score"),
                  colClasses = c("character", "character", "numeric"))
  all_ids <- unique(c(sim_dt$i, sim_dt$j))
  sim_dt  <- rbindlist(list(
    sim_dt,
    data.table(i = sim_dt$j, j = sim_dt$i, score = sim_dt$score),
    data.table(i = all_ids,  j = all_ids,  score = 1.0)
  ))
  setorder(sim_dt, -score)
  sim_dt <- unique(sim_dt, by = c("i", "j"))
  sim_dt
}

# ── Compute pairwise similarity with vsearch --allpairs_global ────────────────
# Output is symmetrised and self-similarities are added before returning.
# ignore.stdout / ignore.stderr suppress per-worker vsearch logging noise.
# tmp_dir is used so workers writing concurrent subset FASTAs do not collide
# with any files in the main output directory.

compute_sim <- function(fasta_file, n_threads, tmp_dir) {
  if (!dir.exists(tmp_dir)) dir.create(tmp_dir, recursive = TRUE)
  vsearch_out <- file.path(tmp_dir, paste0(basename(fasta_file), ".vsearch.txt"))

  vsearch_cmd <- sprintf(
    "vsearch --allpairs_global '%s' --acceptall --userout '%s' --userfields query+target+id --threads %d",
    fasta_file, vsearch_out, n_threads
  )
  if (system(vsearch_cmd, ignore.stdout = TRUE, ignore.stderr = TRUE) != 0) {
    stop("vsearch failed for: ", fasta_file)
  }
  if (!file.exists(vsearch_out) || file.info(vsearch_out)$size == 0) {
    stop("vsearch produced no output for: ", fasta_file)
  }

  # Parse: query <tab> target <tab> percent_identity
  b <- fread(vsearch_out, sep = "\t", header = FALSE,
             col.names = c("qid", "sid", "pident"),
             colClasses = c("character", "character", "numeric"))
  b[, score := round(pident / 100, 4)]

  all_ids <- unique(c(b$qid, b$sid))
  sim_dt  <- rbindlist(list(
    data.table(i = b$qid,   j = b$sid,   score = b$score),
    data.table(i = b$sid,   j = b$qid,   score = b$score),
    data.table(i = all_ids, j = all_ids, score = 1.0)
  ))
  setorder(sim_dt, -score)
  sim_dt <- unique(sim_dt, by = c("i", "j"))

  unlink(vsearch_out)
  sim_dt
}

# ── Save similarity matrix ────────────────────────────────────────────────────

save_sim <- function(sim_dt, path) {
  fwrite(sim_dt, path, sep = " ", col.names = FALSE, quote = FALSE)
}

# ── Write a subset of sequences to a temporary FASTA file ─────────────────────
# Fixed version: uses h_idx[p] to look up the line number of the p-th header,
# avoiding the off-by-one bug present in the original predict.R where
# `h` was incorrectly set to `p` rather than `h_idx[p]`.

write_subset_fasta <- function(fasta_file, ids, out_path) {
  lines    <- readLines(fasta_file)
  h_idx    <- which(startsWith(lines, ">"))
  h_ids    <- sub("^>([^ ]+).*", "\\1", lines[h_idx])
  keep_pos <- which(h_ids %in% ids)
  keep_lines <- unlist(lapply(seq_along(keep_pos), function(k) {
    p   <- keep_pos[k]
    h   <- h_idx[p]
    end <- if (p < length(h_idx)) h_idx[p + 1L] - 1L else length(lines)
    h:end
  }))
  writeLines(lines[keep_lines], out_path)
}

# ── Load the classification table ─────────────────────────────────────────────

load_classification <- function(class_file, ranks, higher_ranks, id_col = "id") {
  cat("[predict] Loading classification from:", class_file, "\n")
  cls <- fread(class_file, sep = "\t", header = TRUE, quote = "",
               fill = TRUE, na.strings = "", check.names = FALSE,
               data.table = FALSE)
  required <- unique(c(id_col, ranks, higher_ranks))
  missing  <- setdiff(required, colnames(cls))
  if (length(missing) > 0) {
    stop("Missing columns in classification file: ", paste(missing, collapse = ", "))
  }
  cls
}

# ── Remove sequences with unidentified / uncertain taxonomy ───────────────────
# Delegates to is_identified() from utils.R for logic centralised filtering.

filter_unidentified_rows <- function(cls, rank) {
  cls[is_identified(cls[[rank]]), ]
}

# ── Generate datasets (global or local by higher-rank taxon) ──────────────────
# Returns a named list of character vectors of sequence IDs.
# Global mode: one dataset named "All" containing all valid sequences.
# Local mode: one dataset per unique value of each higher_rank column.

generate_datasets <- function(seq_ids, cls, rank, higher_ranks,
                              id_col = "id", max_seq_no = 20000) {
  cls_valid <- filter_unidentified_rows(cls[cls[[id_col]] %in% seq_ids, ], rank)

  if (length(higher_ranks) == 0) {
    # Global: one dataset with all valid sequences
    ids <- cls_valid[[id_col]]
    if (max_seq_no > 0 && length(ids) > max_seq_no) ids <- sample(ids, max_seq_no)
    return(list(All = ids))
  }

  # Local: one dataset per higher-rank taxon group
  datasets <- list()
  for (hr in higher_ranks) {
    cls_hr <- cls_valid[
      !is.na(cls_valid[[hr]]) &
      cls_valid[[hr]] != "" &
      cls_valid[[hr]] != "unidentified",
    ]
    groups <- split(cls_hr[[id_col]], cls_hr[[hr]])
    for (grp_name in names(groups)) {
      grp_ids <- groups[[grp_name]]
      if (max_seq_no > 0 && length(grp_ids) > max_seq_no) {
        grp_ids <- sample(grp_ids, max_seq_no)
      }
      datasets[[grp_name]] <- grp_ids
    }
  }
  datasets
}

# ── Group sequences by taxonomic class ────────────────────────────────────────
# Returns a named list of character vectors (class name -> member IDs).

load_classes <- function(seq_ids, cls, rank, id_col = "id") {
  cls_sub <- filter_unidentified_rows(cls[cls[[id_col]] %in% seq_ids, ], rank)
  split(cls_sub[[id_col]], cls_sub[[rank]])
}

# ── Proportion of the dominant (largest) class ────────────────────────────────

max_proportion <- function(classes) {
  n <- sum(lengths(classes))
  if (n == 0) return(1)
  round(max(lengths(classes)) / n, 4)
}

# ── Predict optimal threshold for one dataset ─────────────────────────────────
# Uses a **descending threshold sweep with union-find** for efficiency.
#
# Old approach (O(T × E) total):
#   For each of T threshold steps, re-filter the full edge list, rebuild an
#   adjacency list, run BFS, and compute F-measure with nested intersect() loops.
#   Most of the work is redundant because consecutive thresholds share >99% of
#   their edge sets.
#
# New approach (O(E log E + T × C × K) total, where C = classes, K = clusters):
#   1. Sort edges by score descending (once).
#   2. Sweep thresholds from high -> low.  At each step, add only the *new*
#      edges whose score falls into [t, t+step) and merge components via
#      union-find with path compression + union by rank -> O(α(N)) per edge.
#   3. Maintain a contingency table (cluster root × class) incrementally.
#      On merge, just add the two root vectors -> O(C) per merge.
#   4. Compute F-measure from the contingency table -> O(active_roots × C) per
#      threshold, much cheaper than set intersections.
#
# The result is identical to the old approach but typically 10–30× faster for
# datasets with many threshold steps, because the expensive per-threshold work
# (edge filtering, BFS, intersect) is eliminated entirely.
#
# Caches F-measures from a previous run (in `existing`) to allow incremental
# computation when adding new threshold ranges.
# When the best F-measure is achieved by multiple thresholds, the middle
# threshold of the tied range is selected (rather than the lowest).
# verbose = FALSE suppresses per-threshold output; used inside parallel workers.

predict_dataset <- function(dataset_name, seq_ids, classes, sim_dt,
                            start_t, end_t, step_t,
                            existing = list(), redo = FALSE,
                            verbose = TRUE) {
  # Restore any previously computed F-measures for this dataset
  saved_fm <- if ("fmeasures"  %in% names(existing)) existing$fmeasures   else list()

  # Pre-filter the similarity matrix to only this dataset's sequences
  sub_sim <- sim_dt[i %chin% seq_ids & j %chin% seq_ids]

  if (nrow(sub_sim) == 0) {
    if (verbose) cat("[predict] No similarity data available for:", dataset_name, "\n")
    return(list(error = TRUE))
  }

  # ── Map string IDs -> integers for fast union-find ──────────────────────────
  n         <- length(seq_ids)
  id_to_int <- setNames(seq_len(n), seq_ids)

  # ── Class membership vectors ────────────────────────────────────────────────
  n_classes  <- length(classes)
  class_of   <- integer(n)       # class_of[i] = class index (0 if unclassified)
  class_size <- integer(n_classes)
  for (ci in seq_along(classes)) {
    idx <- id_to_int[classes[[ci]]]
    idx <- idx[!is.na(idx)]
    class_of[idx] <- ci
    class_size[ci] <- length(idx)
  }
  total_n <- sum(class_size)

  # ── Extract unique undirected edges, sort descending by score ───────────────
  edges_dt <- sub_sim[i != j]
  edges_dt[, `:=`(ii = id_to_int[i], jj = id_to_int[j])]
  # Canonicalise direction so each pair appears once (ii < jj)
  swap <- edges_dt$ii > edges_dt$jj
  if (any(swap)) {
    tmp_ii <- edges_dt$ii[swap]
    edges_dt[swap, `:=`(ii = jj, jj = tmp_ii)]
  }
  edges_dt <- unique(edges_dt, by = c("ii", "jj"))
  setorder(edges_dt, -score)

  e_ii    <- edges_dt$ii
  e_jj    <- edges_dt$jj
  e_score <- edges_dt$score
  n_edges <- length(e_score)

  # ── Union-Find arrays ──────────────────────────────────────────────────────
  uf_parent <- seq_len(n)
  uf_rnk    <- integer(n)
  c_size    <- rep(1L, n)     # component size at each root
  is_root   <- rep(TRUE, n)   # track active roots

  # Contingency: root_cc[[r]][c] = count of class-c members in component r
  root_cc <- vector("list", n)
  for (i in seq_len(n)) {
    v <- integer(n_classes)
    if (class_of[i] > 0L) v[class_of[i]] <- 1L
    root_cc[[i]] <- v
  }

  # ── Find with path halving ─────────────────────────────────────────────────
  uf_find <- function(x) {
    while (uf_parent[x] != x) {
      uf_parent[x] <<- uf_parent[uf_parent[x]]
      x <- uf_parent[x]
    }
    x
  }

  # ── Union: merge rb into ra, update contingency ────────────────────────────
  uf_union <- function(a, b) {
    ra <- uf_find(a); rb <- uf_find(b)
    if (ra == rb) return()
    if (uf_rnk[ra] < uf_rnk[rb]) { tmp <- ra; ra <- rb; rb <- tmp }
    uf_parent[rb] <<- ra
    if (uf_rnk[ra] == uf_rnk[rb]) uf_rnk[ra] <<- uf_rnk[ra] + 1L
    c_size[ra]    <<- c_size[ra] + c_size[rb]
    root_cc[[ra]] <<- root_cc[[ra]] + root_cc[[rb]]
    is_root[rb]   <<- FALSE
  }

  # ── Compute F-measure from contingency table ───────────────────────────────
  # For each class c, find the cluster root r maximising
  #   Dice = 2 × count(r,c) / (class_size(c) + component_size(r))
  # Weighted average by class size.
  compute_fm <- function() {
    active <- which(is_root)
    total_f <- 0
    for (ci in seq_len(n_classes)) {
      if (class_size[ci] == 0L) next
      best_dice <- 0
      for (r in active) {
        cnt <- root_cc[[r]][ci]
        if (cnt > 0L) {
          dice <- 2 * cnt / (class_size[ci] + c_size[r])
          if (dice > best_dice) best_dice <- dice
        }
      }
      total_f <- total_f + class_size[ci] * best_dice
    }
    if (total_n == 0L) return(0)
    round(total_f / total_n, 4)
  }

  # ── Sweep thresholds high -> low (add edges as threshold decreases) ─────────
  thresholds_desc <- rev(seq(round(start_t, 4), round(end_t, 4), by = step_t))
  edge_ptr <- 1L   # next edge to process (edges sorted descending by score)
  all_fm   <- list()

  for (t in thresholds_desc) {
    t_str <- sprintf("%.4f", t)

    # Add edges with score >= t that haven't been added yet
    while (edge_ptr <= n_edges && e_score[edge_ptr] >= t - 1e-9) {
      uf_union(e_ii[edge_ptr], e_jj[edge_ptr])
      edge_ptr <- edge_ptr + 1L
    }

    if (!redo && t_str %in% names(saved_fm)) {
      all_fm[[t_str]] <- saved_fm[[t_str]]
    } else {
      fmeasure <- compute_fm()
      all_fm[[t_str]] <- fmeasure
      saved_fm[[t_str]] <- fmeasure
    }
  }

  # ── Collect results in ascending threshold order (matches original output) ─
  thresholds_asc <- seq(round(start_t, 4), round(end_t, 4), by = step_t)
  thresholds <- numeric(length(thresholds_asc))
  fmeasures  <- numeric(length(thresholds_asc))

  for (k in seq_along(thresholds_asc)) {
    t     <- thresholds_asc[k]
    t_str <- sprintf("%.4f", t)
    thresholds[k] <- t
    fmeasures[k]  <- all_fm[[t_str]] %||% 0
  }

  # Select the middle threshold among all thresholds tied at the best F-measure
  best_f      <- max(fmeasures)
  tied_idx    <- which(fmeasures == best_f)
  mid_pos     <- tied_idx[ceiling(length(tied_idx) / 2)]
  opt_t       <- thresholds[mid_pos]

  if (verbose) {
    if (length(tied_idx) > 1) {
      cat(sprintf("[predict] %s: F=%.4f tied over %d thresholds (%.4f\u2013%.4f), selecting middle: %.4f\n",
                  dataset_name, best_f, length(tied_idx),
                  thresholds[tied_idx[1]], thresholds[tied_idx[length(tied_idx)]], opt_t))
    } else {
      cat(sprintf("[predict] %s: optimal cutoff=%.4f  F-measure=%.4f\n",
                  dataset_name, opt_t, best_f))
    }
  }

  list(
    error          = FALSE,
    opt_t          = opt_t,
    best_f         = best_f,
    thresholds     = thresholds,
    fmeasures      = fmeasures,
    fmeasures_dict = saved_fm
  )
}

# ── Process one dataset (called inside each future_map worker) ────────────────
# Encapsulates all per-dataset work: class loading, max-proportion check,
# optional per-dataset vsearch run, and the threshold sweep.
# vsearch_threads: threads given to vsearch within this worker; set to
#   floor(n_cpus / effective_workers) to avoid CPU oversubscription.
# tmp_fasta uses tmp_dir (not output_dir) so parallel workers do not collide
#   with each other or with files in the user-visible output directory.
#
# vsearch_threads is determined per-dataset: datasets with >= large_threshold
# sequences get multi-threaded vsearch; smaller datasets use 1 thread so we
# can run more workers concurrently and keep all CPUs busy during the
# (single-threaded) R threshold sweep.

process_dataset <- function(dataset_name, ds_ids, cls_df, rank, id_col,
                            sim_dt, fasta_file, output_dir,
                            start_t, end_t, step_t, redo, existing,
                            max_prop_limit, vsearch_threads, tmp_dir) {
  classes  <- load_classes(ds_ids, cls_df, rank, id_col)
  n_seqs   <- length(ds_ids)
  n_groups <- length(classes)
  max_prop <- max_proportion(classes)

  if (n_groups < 2) {
    return(list(
      skip = TRUE,
      msg  = sprintf("[predict] Skipping '%s': only %d group(s).",
                     dataset_name, n_groups)
    ))
  }
  if (max_prop >= max_prop_limit && max_prop_limit < 1.0) {
    return(list(
      skip = TRUE,
      msg  = sprintf("[predict] Skipping '%s': dominant group proportion %.4f >= %.4f.",
                     dataset_name, max_prop, max_prop_limit)
    ))
  }

  # Use the pre-loaded global sim_dt when available;
  # otherwise compute similarity on the fly for this subset via vsearch.
  working_sim <- sim_dt
  if (nrow(working_sim) == 0) {
    tmp_fasta <- file.path(tmp_dir,
                           paste0(gsub("[^A-Za-z0-9_]", "_", dataset_name), "_subset.fasta"))
    write_subset_fasta(fasta_file, ds_ids, tmp_fasta)
    working_sim <- compute_sim(tmp_fasta, vsearch_threads, tmp_dir)
    unlink(tmp_fasta)
  }

  result <- predict_dataset(
    dataset_name, ds_ids, classes, working_sim,
    start_t, end_t, step_t, existing, redo,
    verbose = FALSE
  )

  list(
    skip         = FALSE,
    dataset_name = dataset_name,
    result       = result,
    n_seqs       = n_seqs,
    n_groups     = n_groups,
    max_prop     = max_prop
  )
}

# ── Save prediction results ───────────────────────────────────────────────────
# Writes three files:
#   output_file  — full JSON with F-measure traces at every threshold
#   cutoffs_file — filtered JSON with only the best cutoff per dataset
#   cutoffs_file.txt — tab-delimited plain-text summary for easy parsing

save_results <- function(prediction_dict, output_file, cutoffs_file,
                         min_group_no, min_seq_no, max_prop_limit, min_cutoff) {
  # Full prediction file (includes per-threshold F-measure traces)
  write(toJSON(prediction_dict, auto_unbox = TRUE, pretty = TRUE), output_file)

  # Cutoffs file: filtered to datasets that pass quality thresholds
  final <- prediction_dict
  for (rank in names(final)) {
    to_rm <- character(0)
    for (dname in names(final[[rank]])) {
      d       <- final[[rank]][[dname]]
      groupno <- d[["group number"]]    %||% 0
      seqno   <- d[["sequence number"]] %||% 0
      maxprop <- d[["max proportion"]]  %||% 0
      cutoff  <- d[["cut-off"]]         %||% 0

      if (groupno < min_group_no || seqno < min_seq_no ||
          maxprop > max_prop_limit || cutoff < min_cutoff) {
        to_rm <- c(to_rm, dname)
      } else {
        # Strip the bulky fmeasures trace from the cutoffs-only output
        final[[rank]][[dname]][["fmeasures"]] <- NULL
      }
    }
    for (dn in to_rm) final[[rank]][[dn]] <- NULL
  }

  write(toJSON(final, auto_unbox = TRUE, pretty = TRUE), cutoffs_file)

  # Tab-delimited plain-text summary
  txt_file <- paste0(cutoffs_file, ".txt")
  header   <- paste(c("rank", "higher_rank", "dataset", "cut-off", "confidence",
                       "sequence number", "group number", "max proportion"),
                    collapse = "\t")
  rows <- header
  for (rank in names(final)) {
    for (dname in names(final[[rank]])) {
      d <- final[[rank]][[dname]]
      rows <- c(rows, paste(
        rank,
        d[["higher_rank"]]     %||% "global",
        dname,
        d[["cut-off"]]         %||% 0,
        d[["confidence"]]      %||% 0,
        d[["sequence number"]] %||% 0,
        d[["group number"]]    %||% 0,
        d[["max proportion"]]  %||% 0,
        sep = "\t"
      ))
    }
  }
  writeLines(rows, txt_file)

  cat("[predict] Full prediction (with F-measure traces):", output_file, "\n")
  cat("[predict] Cutoffs JSON:", cutoffs_file, "\n")
  cat("[predict] Cutoffs text:", txt_file, "\n")
}

# ══════════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════════

rank_list        <- trimws(strsplit(rank_str, ",")[[1]])
higher_rank_list <- if (nchar(higher_rank_str) > 0) trimws(strsplit(higher_rank_str, ",")[[1]]) else character(0)

# ── Load existing prediction to allow incremental / cached runs ───────────────
prediction_dict <- list()
if (file.exists(prediction_file)) {
  cat("[predict] Loading existing prediction from:", prediction_file, "\n")
  prediction_dict <- fromJSON(prediction_file, simplifyVector = FALSE)
}

# ── Load sequence IDs from FASTA ─────────────────────────────────────────────
cat("[predict] Loading sequences from:", fasta_file, "\n")
fasta_info <- load_fasta_ids(fasta_file)
seq_ids    <- fasta_info$ids
cat("[predict]", length(seq_ids), "sequences loaded.\n")

# ── Resolve similarity matrix path ────────────────────────────────────────────
# If --sim was provided use it directly; otherwise default to the output dir.
sim_file_path <- if (nchar(sim_file_arg) > 0) sim_file_arg else
                   file.path(output_dir, paste0(prefix, ".sim"))

# ── Load or compute the similarity matrix ─────────────────────────────────────
sim_dt <- data.table(i = character(0), j = character(0), score = numeric(0))

if (end_t >= start_t && end_t > 0) {
  if (file.exists(sim_file_path)) {
    sim_dt <- load_sim(sim_file_path)
  } else if (length(higher_rank_list) == 0) {
    # Global prediction: compute similarity once for all sequences
    if (length(seq_ids) <= max_seq_no) {
      cat("[predict] No .sim file found. Computing similarity with vsearch...\n")
      sim_dt <- compute_sim(fasta_file, n_cpus, tmp_dir)
      save_sim(sim_dt, sim_file_path)
      cat("[predict] Similarity matrix saved to:", sim_file_path, "\n")
    } else {
      cat("[predict] Dataset too large to auto-compute similarity.",
          "Run sim_vsearch.R first and supply --sim.\n")
    }
  } else {
    cat("[predict] No .sim file found. For local prediction, compute similarity first",
        "and supply --sim.\n",
        "Alternatively, similarity will be computed per-dataset via vsearch.\n")
  }
}

# ── Load classification table ─────────────────────────────────────────────────
cls_df <- load_classification(class_file, rank_list, higher_rank_list, id_col)

# ══════════════════════════════════════════════════════════════════════════════
# Rank loop
# ══════════════════════════════════════════════════════════════════════════════

for (rank in rank_list) {
  cat(sprintf("\n[predict] ══════ Rank: %s ══════\n", rank))

  if (!rank %in% names(prediction_dict)) prediction_dict[[rank]] <- list()
  pred_datasets <- prediction_dict[[rank]]

  datasets <- generate_datasets(seq_ids, cls_df, rank, higher_rank_list,
                                id_col, max_seq_no)

  if (length(datasets) == 0) {
    cat("[predict] No valid datasets found for rank:", rank, "\n")
    next
  }
  if (end_t < start_t || end_t == 0) {
    cat("[predict] Threshold range not specified — skipping computation.\n")
    next
  }

  n_datasets  <- length(datasets)
  max_ds_size <- max(sapply(datasets, length))

  # ── Adaptive threading (two-batch strategy) ──────────────────────────────
  # Problem:  The old approach picked a single worker×thread split based on
  #   the *largest* dataset.  When one dataset is large but most are tiny,
  #   this caps workers at n_cpus/8 ≈ 10 — but the R threshold sweep inside
  #   each worker is single-threaded, so only ~10 CPUs stay busy.
  # Solution: Split datasets into "large" (>= large_threshold seqs, benefits
  #   from multi-threaded vsearch) and "small" (< large_threshold, vsearch
  #   finishes instantly so maximise R-worker concurrency).
  #   Batch 1 — large datasets:  fewer workers × more vsearch threads.
  #   Batch 2 — small datasets:  many workers × 1 vsearch thread.
  #   When a pre-computed sim file was loaded vsearch is never called, so all
  #   datasets go into the small (high-worker) batch.
  # ── Scaled threading: vsearch threads scale linearly with dataset size ──
  #
  # Rationale: large datasets need more vsearch threads for the O(N²) pairwise
  # alignment but that means fewer concurrent workers, keeping total CPU usage
  # ≈ n_cpus while reducing per-worker memory (the OOM-kill root cause).
  #
  # Scaling schedule (linear interpolation):
  #   < 400 seqs  -> 1  vsearch thread  (tiny; maximise workers)
  #   400 seqs    -> 2  vsearch threads
  #   ≈12 500 seqs -> n_cpus vsearch threads (single worker gets all CPUs)
  #
  # The schedule rounds to the nearest value in {1,2,4,8,16,24,32,40,...}
  # so workers divide evenly into n_cpus.

  has_sim  <- nrow(sim_dt) > 0
  ds_sizes <- sapply(datasets, length)

  # Permitted thread counts (powers-of-2 only, up to n_cpus). Every tier
  # evenly divides n_cpus so all CPUs stay busy (e.g. 32/4 = 8 workers).
  thread_steps <- 2L^(0:floor(log2(n_cpus)))
  thread_steps <- unique(sort(c(thread_steps, n_cpus)))
  thread_steps <- thread_steps[thread_steps <= n_cpus]

  # Linear interpolation helper
  SCALE_MIN_SIZE    <- 400
  SCALE_MAX_SIZE    <- 12500
  SCALE_MIN_THREADS <- 2L
  SCALE_MAX_THREADS <- n_cpus

  compute_threads <- function(ds_size) {
    if (ds_size < SCALE_MIN_SIZE) return(1L)
    raw <- SCALE_MIN_THREADS +
      (ds_size - SCALE_MIN_SIZE) *
      (SCALE_MAX_THREADS - SCALE_MIN_THREADS) /
      (SCALE_MAX_SIZE - SCALE_MIN_SIZE)
    raw <- max(SCALE_MIN_THREADS, min(SCALE_MAX_THREADS, raw))
    # Snap to the nearest permitted thread step
    as.integer(thread_steps[which.min(abs(thread_steps - raw))])
  }

  # Assign threads per dataset (1 thread if pre-computed sim available)
  if (has_sim) {
    thread_assignments <- setNames(rep(1L, length(ds_sizes)), names(ds_sizes))
  } else {
    thread_assignments <- setNames(
      vapply(ds_sizes, compute_threads, integer(1)),
      names(ds_sizes)
    )
  }

  # Group datasets into tiers by their thread count
  tier_levels <- sort(unique(thread_assignments))

  cat(sprintf("[predict] %d dataset(s), largest=%d seqs, %d tier(s) (threads scaled %d->%d).\n",
              n_datasets, max_ds_size, length(tier_levels),
              min(tier_levels), max(tier_levels)))
  for (thr in tier_levels) {
    tier_ns <- names(thread_assignments[thread_assignments == thr])
    tier_workers <- max(1L, floor(n_cpus / thr))
    cat(sprintf(
      "[predict]   Tier %2d-thread: %3d dataset(s) -> %d worker(s) × %d vsearch thread(s).\n",
      thr, length(tier_ns), min(length(tier_ns), tier_workers), thr
    ))
  }

  # ── Helper: build argument list for a set of datasets ─────────────────────
  # Sort largest-first so big jobs start immediately; smaller tasks fill in
  # gaps as workers finish, reducing tail latency.
  make_args <- function(ds_names) {
    ds_names <- ds_names[order(ds_sizes[ds_names], decreasing = TRUE)]
    lapply(ds_names, function(dn) {
      list(
        dataset_name = dn,
        ds_ids       = datasets[[dn]],
        existing     = if (dn %in% names(pred_datasets)) pred_datasets[[dn]] else list()
      )
    })
  }

  # ── Helper: run one batch via future_map (with sequential fallback) ──────
  #
  # When parallel execution crashes (typically because the Linux OOM killer
  # terminates a forked child), the batch is automatically re-run
  # sequentially so that (a) the specific failing datasets are identified
  # by name and (b) the error message / reason is captured in the log.

  run_batch <- function(args_list, batch_workers, batch_vsearch_threads) {
    if (length(args_list) == 0) return(list())

    # ── Inner helper: process one dataset with error trapping ─────────────
    do_one <- function(a, vthreads) {
      tryCatch(
        process_dataset(
          dataset_name    = a$dataset_name,
          ds_ids          = a$ds_ids,
          cls_df          = cls_df,
          rank            = rank,
          id_col          = id_col,
          sim_dt          = sim_dt,
          fasta_file      = fasta_file,
          output_dir      = output_dir,
          start_t         = start_t,
          end_t           = end_t,
          step_t          = step_t,
          redo            = redo,
          existing        = a$existing,
          max_prop_limit  = max_prop_limit,
          vsearch_threads = vthreads,
          tmp_dir         = tmp_dir
        ),
        error = function(e) {
          list(
            skip         = FALSE,
            dataset_name = a$dataset_name,
            result       = list(error    = TRUE,
                                error_msg = conditionMessage(e)),
            n_seqs       = length(a$ds_ids),
            n_groups     = NA_integer_,
            max_prop     = NA_real_
          )
        }
      )
    }

    # ── Attempt parallel execution ───────────────────────────────────────
    if (run_parallel) {
      if (.Platform$OS.type == "unix") {
        plan(multicore,    workers = min(length(args_list), batch_workers))
      } else {
        plan(multisession, workers = min(length(args_list), batch_workers))
      }

      parallel_ok <- tryCatch({
        result <- future_map(
          args_list,
          function(a) do_one(a, batch_vsearch_threads),
          .options  = furrr_options(seed = NULL, chunk_size = 1L),
          .progress = FALSE
        )
        TRUE
      }, error = function(e) {
        cat(sprintf(
          "\n[predict] WARNING: Parallel batch crashed — %s\n",
          conditionMessage(e)
        ))
        cat("[predict]   This is typically caused by the OS OOM (out-of-memory) killer\n")
        cat("[predict]   terminating forked child processes.\n")
        cat("[predict]   Falling back to sequential processing to identify failing datasets...\n\n")
        FALSE
      })

      if (parallel_ok) return(result)
    }

    # ── Sequential fallback (or non-parallel mode) ───────────────────────
    plan(sequential)
    lapply(args_list, function(a) {
      cat(sprintf("[predict]   Sequential: processing '%s' (%d seqs)...\n",
                  a$dataset_name, length(a$ds_ids)))
      res <- do_one(a, n_cpus)   # give all threads to the single dataset
      if (!is.null(res$result$error_msg)) {
        cat(sprintf("[predict]   >>> FAILED '%s': %s\n",
                    a$dataset_name, res$result$error_msg))
      }
      res
    })
  }

  # ── Run tiered batches (smallest datasets first -> largest last) ──────────
  results_list <- list()
  for (thr in tier_levels) {
    tier_names   <- names(thread_assignments[thread_assignments == thr])
    tier_workers <- max(1L, floor(n_cpus / thr))
    tier_results <- run_batch(make_args(tier_names), tier_workers, thr)
    results_list <- c(results_list, tier_results)
  }

  # ── Collect results and update pred_datasets ──────────────────────────────
  # For small runs (≤ 100 datasets) print one line per result so test output
  # is readable.  For large runs print a progress summary every ~5% to avoid
  # flooding the SLURM log with tens of thousands of lines.
  verbose_per_dataset <- n_datasets <= 100L
  report_every        <- if (verbose_per_dataset) 1L else max(1L, floor(n_datasets / 20L))
  n_done    <- 0L
  n_cutoffs <- 0L
  n_skipped <- 0L
  n_errors  <- 0L

  for (res in results_list) {
    n_done <- n_done + 1L
    if (res$skip) {
      n_skipped <- n_skipped + 1L
      if (verbose_per_dataset) cat(res$msg, "\n")
    } else {
      dn     <- res$dataset_name
      result <- res$result
      if (!result$error) {
        n_cutoffs <- n_cutoffs + 1L
        if (verbose_per_dataset) {
          cat(sprintf("[predict] %s: cutoff=%.4f  F=%.4f  (seqs=%d  groups=%d)\n",
                      dn, result$opt_t, result$best_f, res$n_seqs, res$n_groups))
        }
        pred_datasets[[dn]] <- list(
          "higher_rank"     = if (nchar(higher_rank_str) > 0) higher_rank_str else "global",
          "cut-off"         = result$opt_t,
          "confidence"      = result$best_f,
          "sequence number" = res$n_seqs,
          "group number"    = res$n_groups,
          "max proportion"  = res$max_prop,
          "fmeasures"       = result$fmeasures_dict
        )
      } else {
        n_errors <- n_errors + 1L
        err_detail <- if (!is.null(result$error_msg)) result$error_msg
                      else "no similarity data or unknown error"
        cat(sprintf("[predict] ERROR: dataset '%s' failed \u2014 %s\n", dn, err_detail))
      }
    }
    if (!verbose_per_dataset && (n_done %% report_every == 0L || n_done == n_datasets)) {
      cat(sprintf("[predict]   %d/%d datasets  cutoffs: %d  skipped: %d\n",
                  n_done, n_datasets, n_cutoffs, n_skipped))
    }
  }
  cat(sprintf("[predict] Rank %s: %d cutoff(s) from %d dataset(s) (%d skipped, %d error(s)).\n",
              rank, n_cutoffs, n_datasets, n_skipped, n_errors))

  # ── Log a summary of all failed datasets (if any) ──────────────────────
  if (n_errors > 0L) {
    failed_names <- character(0)
    for (res in results_list) {
      if (!isTRUE(res$skip) && isTRUE(res$result$error)) {
        failed_names <- c(failed_names, res$dataset_name)
      }
    }
    cat(sprintf("[predict] *** %d dataset(s) FAILED for rank '%s':\n", n_errors, rank))
    for (fn in failed_names) cat(sprintf("[predict]     - %s\n", fn))
    cat("[predict] *** Likely cause: OS OOM killer terminated the worker process.\n")
    cat("[predict]     Consider reducing --max_seq_no, requesting more memory in SLURM,\n")
    cat("[predict]     or processing these datasets with --run_parallel no.\n")
  }

  prediction_dict[[rank]] <- pred_datasets
}

# ── Save ──────────────────────────────────────────────────────────────────────

has_results <- any(sapply(prediction_dict, function(r) length(r) > 0))

if (has_results) {
  cutoffs_file <- sub("\\.predicted$", ".cutoffs.json", prediction_file)
  save_results(prediction_dict, prediction_file, cutoffs_file,
               min_group_no, min_seq_no, max_prop_limit, min_cutoff)
  cat(sprintf(
    "\n[predict] Done. Cutoffs reported for datasets with >= %d sequences,",
    min_seq_no
  ))
  cat(sprintf(
    " >= %d groups, max proportion < %.2f, and cutoff >= %.4f.\n",
    min_group_no, max_prop_limit, min_cutoff
  ))
} else {
  cat("[predict] No prediction results to save.\n")
}

plan(sequential)  # restore default plan on exit
