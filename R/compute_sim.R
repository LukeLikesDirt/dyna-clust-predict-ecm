#!/usr/bin/env Rscript
# sim_vsearch.R — Pre-compute pairwise vsearch similarity matrix
#
# Usage:
#   Rscript sim_vsearch.R --input sequences.fasta [options]
#
# Required packages: optparse
# Requires vsearch on PATH: https://github.com/torognes/vsearch
#
# Output:
#   <out>/<input-basename>.sim   — space-delimited: id1 id2 score (no header)
#
# Note: This script must be run from the project root directory.

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))

# ── Arguments ─────────────────────────────────────────────────────────────────

option_list <- list(
  make_option(c("-i", "--input"),
              type = "character", metavar = "FILE",
              help = "Input FASTA file [required]"),
  make_option(c("-o", "--out"),
              type = "character", default = ".", metavar = "DIR",
              help = "Output directory [default: %default]"),
  make_option("--min_sim",
              type = "double", default = 0, metavar = "NUM",
              help = "Minimum similarity score to write to output [default: %default]"),
  make_option("--iddef",
              type = "integer", default = 2L, metavar = "0-4",
              help = paste("vsearch --iddef pairwise identity definition.",
                           "Default 2 (terminal gaps excluded), as in R/predict.R.",
                           "[default: %default]")),
  make_option("--min_overlap",
              type = "integer", default = 0L, metavar = "BP",
              help = paste("Drop pairs whose alignment, terminal gaps excluded, is shorter",
                           "than this many bp (vsearch --mincols). See R/predict.R",
                           "--min_overlap. Default 0 (off). [default: %default]")),
  make_option("--n_cpus",
              type = "integer", default = max(1L, parallel::detectCores() - 1L), metavar = "INT",
              help = "CPU threads for vsearch [default: all cores minus one]"),
  make_option("--tmp_dir",
              type = "character", default = "./tmp", metavar = "DIR",
              help = "Directory for temporary vsearch output [default: %default]")
)

opt <- parse_args(
  OptionParser(
    option_list  = option_list,
    usage        = "%prog [options] -i sequences.fasta",
    description  = paste(
      "Compute a pairwise vsearch similarity matrix using --allpairs_global.",
      "Similarity = pident / 100 (global alignment identity, --iddef 2 by default:",
      "terminal gaps excluded; --min_overlap drops pairs that share too few bases).",
      "The matrix is symmetrised; the maximum score is kept for each pair."
    )
  )
)

if (is.null(opt$input))      stop("--input is required. Run with --help for usage.")
if (!file.exists(opt$input)) stop("FASTA file not found: ", opt$input)
if (is.na(opt$iddef) || !opt$iddef %in% 0:4) stop("--iddef must be 0, 1, 2, 3 or 4.")
if (is.na(opt$min_overlap) || opt$min_overlap < 0) {
  stop("--min_overlap must be a non-negative number of bp (0 disables it).")
}

fasta_file  <- opt$input
output_dir  <- opt$out
min_sim     <- opt$min_sim
iddef       <- opt$iddef
min_overlap <- opt$min_overlap
n_cpus      <- opt$n_cpus
tmp_dir    <- opt$tmp_dir

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
if (!dir.exists(tmp_dir))    dir.create(tmp_dir,    recursive = TRUE)

# ── Helpers ───────────────────────────────────────────────────────────────────

strip_ext <- function(path) sub("\\.[^./]+$", "", basename(path))

# ── vsearch all-vs-all ────────────────────────────────────────────────────────

compute_vsearch_sim <- function(fasta_file, n_cpus, tmp_dir, iddef = 2L, min_overlap = 0L) {
  vsearch_out <- file.path(tmp_dir, paste0(basename(fasta_file), ".vsearch.txt"))

  # --acceptall switches off every accept/reject criterion, --mincols included,
  # so a minimum overlap accepts any identity (--id 0) instead. Rejected pairs
  # are absent from vsearch's output, so the node set comes from the FASTA: a
  # sequence whose every pair is rejected stays as a self-only node. See
  # compute_sim() in R/predict.R.
  gated  <- min_overlap > 0
  accept <- if (gated) sprintf("--id 0 --mincols %d", as.integer(min_overlap)) else "--acceptall"
  vsearch_cmd <- sprintf(
    "vsearch --allpairs_global '%s' %s --userout '%s' --userfields query+target+id --iddef %d --threads %d",
    fasta_file, accept, vsearch_out, iddef, n_cpus
  )
  cat("[sim] Running vsearch --allpairs_global...\n")
  if (system(vsearch_cmd) != 0) stop("vsearch failed. Is vsearch installed and on PATH?")

  no_hits <- !file.exists(vsearch_out) || file.info(vsearch_out)$size == 0
  if (no_hits && !gated) {
    stop("vsearch produced no output. Check that the FASTA file contains valid sequences.")
  }

  cat("[sim] Parsing vsearch output...\n")
  b <- if (no_hits) {
    data.table(qid = character(0), sid = character(0), pident = numeric(0))
  } else {
    fread(vsearch_out, sep = "\t", header = FALSE,
          col.names = c("qid", "sid", "pident"),
          nThread = n_cpus)
  }
  b[, score := round(pident / 100, 4)]
  b[, pident := NULL]

  all_ids <- if (gated) {
    lines <- readLines(fasta_file)
    sub("^>([^ ]+).*", "\\1", lines[startsWith(lines, ">")])
  } else {
    unique(c(b$qid, b$sid))
  }
  sim_dt  <- rbindlist(list(
    b[, .(i = qid, j = sid, score)],
    b[, .(i = sid, j = qid, score)],
    data.table(i = all_ids, j = all_ids, score = 1.0)
  ))
  setorder(sim_dt, -score)
  sim_dt <- unique(sim_dt, by = c("i", "j"))

  unlink(vsearch_out)
  sim_dt
}

# ── Main ──────────────────────────────────────────────────────────────────────

cat("[sim] Input:   ", fasta_file, "\n")
cat("[sim] Threads: ", n_cpus, "\n")
cat("[sim] Tmp dir: ", tmp_dir, "\n")
cat("[sim] iddef:   ", iddef, "\n")
cat("[sim] min_overlap: ", min_overlap, "bp\n\n")

sim_df   <- compute_vsearch_sim(fasta_file, n_cpus, tmp_dir, iddef, min_overlap)
out_file <- file.path(output_dir, paste0(strip_ext(fasta_file), ".sim"))

cat("[sim] Saving similarity matrix to:", out_file, "\n")
out_dt <- sim_df[score >= min_sim]
fwrite(out_dt, out_file, sep = " ", quote = FALSE, col.names = FALSE, nThread = n_cpus)

cat(sprintf("[sim] Done. %d entries saved (%.0f unique sequences).\n",
            nrow(out_dt), length(unique(c(out_dt$i, out_dt$j)))))
