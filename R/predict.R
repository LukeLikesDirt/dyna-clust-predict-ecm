#!/usr/bin/env Rscript
# predict.R — Predict optimal similarity cutoffs for DNA barcoding taxonomy.
#
# Unified script merging predict_vsearch.R and predict_vsearch_parallel.R.
# Supports both sequential and parallel (furrr/future) dataset processing,
# controlled by --run_parallel (default: yes).  Parallel/sequential mode does
# NOT affect which packages are loaded — all packages are always loaded.
#
# Uses vsearch --allpairs_global (global alignment) for similarity computation.
# A pre-computed .sim file may be supplied via --sim to skip that step.
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
              help = paste("Higher rank(s) for local prediction, comma-separated (e.g. genus).",
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
              type = "integer", default = 20000L, metavar = "INT",
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
# avoiding the off-by-one bug present in the original predict_vsearch.R where
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
# Uses fread() for fast loading of large TSV files.

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
# Delegates to is_identified() from utils.R so filtering logic is centralised.

filter_unidentified_rows <- function(cls, rank) {
  cls[is_identified(cls[[rank]]), ]
}

# ── Generate datasets (global or local by higher-rank taxon) ──────────────────
# Returns a named list of character vectors of sequence IDs.
# Global mode: one dataset named "All" containing all valid sequences.
# Local mode:  one dataset per unique value of each higher_rank column.

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
# Returns a named list of character vectors (class name → member IDs).

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

# ── Build neighbour list: sequences connected at sim >= threshold ──────────────
# sub_sim must already be filtered to the relevant sequence IDs.
# Uses data.table column-filtering for speed.

build_neighbors <- function(ids, sub_sim, threshold) {
  edges <- sub_sim[score >= threshold, .(i, j)]
  adj   <- split(edges$j, edges$i)
  neighbor_list <- setNames(vector("list", length(ids)), ids)
  for (id in ids) neighbor_list[[id]] <- adj[[id]]
  neighbor_list
}

# ── Connected-components clustering (iterative BFS) ───────────────────────────
# Returns a list of character vectors; each element is the IDs in one cluster.
# Iterative (not recursive) to avoid stack overflow on large datasets.

find_clusters <- function(ids, neighbor_list) {
  visited  <- setNames(logical(length(ids)), ids)
  clusters <- list()

  for (start in ids) {
    if (!visited[start]) {
      queue          <- start
      visited[start] <- TRUE
      members        <- start

      while (length(queue) > 0) {
        current <- queue[1]
        queue   <- queue[-1]

        nbrs <- neighbor_list[[current]]
        nbrs <- nbrs[nbrs %in% ids]
        new  <- nbrs[!visited[nbrs]]

        if (length(new) > 0) {
          visited[new] <- TRUE
          members      <- c(members, new)
          queue        <- c(queue, new)
        }
      }

      clusters[[length(clusters) + 1L]] <- members
    }
  }
  clusters
}

# ── F-measure between true taxonomic classes and predicted clusters ───────────
# Weighted average over classes of the best-matching cluster's Dice coefficient:
#   2 * |intersection| / (|class| + |cluster|)
# Weighted by class size so larger classes contribute more.

compute_fmeasure <- function(classes, clusters) {
  f <- 0
  n <- 0
  for (group in classes) {
    m <- 0
    for (cl in clusters) {
      i_size <- length(intersect(group, cl))
      v      <- 2 * i_size / (length(group) + length(cl))
      if (v > m) m <- v
    }
    n <- n + length(group)
    f <- f + length(group) * m
  }
  if (n == 0) return(0)
  round(f / n, 4)
}

# ── Predict optimal threshold for one dataset ─────────────────────────────────
# Sweeps thresholds from start_t to end_t in steps of step_t.
# Caches F-measures from a previous run (in `existing`) to allow incremental
# computation when adding new threshold ranges.
# verbose = FALSE suppresses per-threshold output; used inside parallel workers
# where many datasets may log concurrently and make output unreadable.

predict_dataset <- function(dataset_name, seq_ids, classes, sim_dt,
                            start_t, end_t, step_t,
                            existing = list(), redo = FALSE,
                            verbose = TRUE) {
  # Restore any previously computed F-measures for this dataset
  saved_fm <- if ("fmeasures"  %in% names(existing)) existing$fmeasures   else list()
  opt_t    <- if ("cut-off"    %in% names(existing)) existing[["cut-off"]] else 0
  best_f   <- if ("confidence" %in% names(existing)) existing$confidence   else 0

  # Reset cached optimum when redo is requested or threshold range has changed
  if (redo || opt_t < start_t || (end_t > 0 && opt_t > end_t)) {
    opt_t  <- 0
    best_f <- 0
  }

  # Pre-filter the similarity matrix to only this dataset's sequences
  # %chin% is data.table's fast character %in%
  sub_sim <- sim_dt[i %chin% seq_ids & j %chin% seq_ids]

  if (nrow(sub_sim) == 0) {
    if (verbose) cat("[predict] No similarity data available for:", dataset_name, "\n")
    return(list(error = TRUE))
  }

  thresholds <- numeric(0)
  fmeasures  <- numeric(0)
  t          <- round(start_t, 4)

  while (t <= end_t + 1e-9) {  # small epsilon handles floating-point rounding
    t_str <- sprintf("%.4f", t)

    if (!redo && t_str %in% names(saved_fm)) {
      fmeasure <- saved_fm[[t_str]]
    } else {
      neighbors <- build_neighbors(seq_ids, sub_sim, t)
      clusters  <- find_clusters(seq_ids, neighbors)
      fmeasure  <- compute_fmeasure(classes, clusters)
      saved_fm[[t_str]] <- fmeasure
    }

    if (fmeasure > best_f || (fmeasure == best_f && t < opt_t)) {
      best_f <- fmeasure
      opt_t  <- t
    }

    thresholds <- c(thresholds, t)
    fmeasures  <- c(fmeasures,  fmeasure)

    if (verbose) cat(sprintf("  threshold=%.4f  F=%.4f\n", t, fmeasure))
    t <- round(t + step_t, 4)
  }

  if (verbose) {
    cat(sprintf("[predict] %s: optimal cutoff=%.4f  F-measure=%.4f\n",
                dataset_name, opt_t, best_f))
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
  header   <- paste(c("Rank", "Dataset", "cut-off", "confidence",
                       "sequence number", "group number", "max proportion"),
                    collapse = "\t")
  rows <- header
  for (rank in names(final)) {
    for (dname in names(final[[rank]])) {
      d <- final[[rank]][[dname]]
      rows <- c(rows, paste(
        rank, dname,
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

  # ── Adaptive threading ────────────────────────────────────────────────────
  # vsearch --allpairs_global is O(N²) in pairwise alignments; for large
  # datasets multi-threading gives near-linear speedup. For small datasets
  # (< 500 seqs, ~25k pairs) single-threaded vsearch is fast enough and
  # maximising the number of parallel dataset workers is more valuable.
  # When a pre-computed sim file was loaded vsearch is not called at all, so
  # maximise workers unconditionally in that case.
  large_threshold <- 500L
  if (nrow(sim_dt) == 0 && max_ds_size >= large_threshold) {
    vsearch_threads   <- min(8L, n_cpus)
    effective_workers <- max(1L, floor(n_cpus / vsearch_threads))
  } else {
    vsearch_threads   <- 1L
    effective_workers <- min(n_datasets, n_cpus)
  }

  cat(sprintf(
    "[predict] %d dataset(s), largest=%d seqs → %d worker(s) × %d vsearch thread(s).\n",
    n_datasets, max_ds_size, min(n_datasets, effective_workers), vsearch_threads
  ))

  # ── Adjust plan to effective_workers for this rank ────────────────────────
  if (run_parallel) {
    if (.Platform$OS.type == "unix") {
      plan(multicore,    workers = min(n_datasets, effective_workers))
    } else {
      plan(multisession, workers = min(n_datasets, effective_workers))
    }
  }

  # ── Build argument list (pass existing results so cached F-measures reused) ─
  args_list <- lapply(names(datasets), function(dn) {
    list(
      dataset_name = dn,
      ds_ids       = datasets[[dn]],
      existing     = if (dn %in% names(pred_datasets)) pred_datasets[[dn]] else list()
    )
  })

  # ── Map over datasets (parallel or sequential via furrr) ──────────────────
  # .progress = FALSE: in non-TTY contexts (SLURM) progressr prints a new line
  # per completed task rather than updating in-place, flooding the log file.
  results_list <- future_map(
    args_list,
    function(a) {
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
        vsearch_threads = vsearch_threads,
        tmp_dir         = tmp_dir
      )
    },
    .options  = furrr_options(seed = NULL),
    .progress = FALSE
  )

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
          "cut-off"         = result$opt_t,
          "confidence"      = result$best_f,
          "sequence number" = res$n_seqs,
          "group number"    = res$n_groups,
          "max proportion"  = res$max_prop,
          "fmeasures"       = result$fmeasures_dict
        )
      } else {
        n_errors <- n_errors + 1L
        cat(sprintf("[predict] WARNING: no similarity data for '%s', skipped.\n", dn))
      }
    }
    if (!verbose_per_dataset && (n_done %% report_every == 0L || n_done == n_datasets)) {
      cat(sprintf("[predict]   %d/%d datasets  cutoffs: %d  skipped: %d\n",
                  n_done, n_datasets, n_cutoffs, n_skipped))
    }
  }
  cat(sprintf("[predict] Rank %s: %d cutoff(s) from %d dataset(s) (%d skipped, %d error(s)).\n",
              rank, n_cutoffs, n_datasets, n_skipped, n_errors))

  prediction_dict[[rank]] <- pred_datasets

  # ── Restore full worker count for the next rank ───────────────────────────
  if (run_parallel) {
    if (.Platform$OS.type == "unix") {
      plan(multicore,    workers = n_cpus)
    } else {
      plan(multisession, workers = n_cpus)
    }
  }
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
