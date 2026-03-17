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
      "Similarity = pident / 100 (global alignment identity, no coverage penalty).",
      "The matrix is symmetrised; the maximum score is kept for each pair."
    )
  )
)

if (is.null(opt$input))      stop("--input is required. Run with --help for usage.")
if (!file.exists(opt$input)) stop("FASTA file not found: ", opt$input)

fasta_file <- opt$input
output_dir <- opt$out
min_sim    <- opt$min_sim
n_cpus     <- opt$n_cpus
tmp_dir    <- opt$tmp_dir

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
if (!dir.exists(tmp_dir))    dir.create(tmp_dir,    recursive = TRUE)

# ── Helpers ───────────────────────────────────────────────────────────────────

strip_ext <- function(path) sub("\\.[^./]+$", "", basename(path))

# ── vsearch all-vs-all ────────────────────────────────────────────────────────

compute_vsearch_sim <- function(fasta_file, n_cpus, tmp_dir) {
  vsearch_out <- file.path(tmp_dir, paste0(basename(fasta_file), ".vsearch.txt"))

  vsearch_cmd <- sprintf(
    "vsearch --allpairs_global '%s' --acceptall --userout '%s' --userfields query+target+id --threads %d",
    fasta_file, vsearch_out, n_cpus
  )
  cat("[sim] Running vsearch --allpairs_global...\n")
  if (system(vsearch_cmd) != 0) stop("vsearch failed. Is vsearch installed and on PATH?")

  if (!file.exists(vsearch_out) || file.info(vsearch_out)$size == 0) {
    stop("vsearch produced no output. Check that the FASTA file contains valid sequences.")
  }

  cat("[sim] Parsing vsearch output...\n")
  b <- fread(vsearch_out, sep = "\t", header = FALSE,
             col.names = c("qid", "sid", "pident"), nThread = n_cpus)
  b[, score := round(pident / 100, 4)]
  b[, pident := NULL]

  all_ids <- unique(c(b$qid, b$sid))
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
cat("[sim] Tmp dir: ", tmp_dir, "\n\n")

sim_df   <- compute_vsearch_sim(fasta_file, n_cpus, tmp_dir)
out_file <- file.path(output_dir, paste0(strip_ext(fasta_file), ".sim"))

cat("[sim] Saving similarity matrix to:", out_file, "\n")
out_dt <- sim_df[score >= min_sim]
fwrite(out_dt, out_file, sep = " ", quote = FALSE, col.names = FALSE, nThread = n_cpus)

cat(sprintf("[sim] Done. %d entries saved (%.0f unique sequences).\n",
            nrow(out_dt), length(unique(c(out_dt$i, out_dt$j)))))
