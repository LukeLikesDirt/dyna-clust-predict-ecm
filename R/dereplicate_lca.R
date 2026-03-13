#!/usr/bin/env Rscript
# dereplicate_lca.R — Dereplicate FASTA and resolve taxonomy with LCA
#
# Groups sequences by exact string identity, keeps the first representative
# per group, and resolves taxonomy to the lowest common ancestor (LCA) by
# walking ranks top-down and stopping at the first rank with conflicting values.
#
# Usage:
#   Rscript dereplicate_lca.R \
#     --fasta_in              data/full_ITS/eukaryome_ITS_filtered.fasta \
#     --fasta_out             data/full_ITS/eukaryome_ITS_dereplicated.fasta \
#     --classification_in     data/full_ITS/eukaryome_ITS_filtered.classification \
#     --classification_out    data/full_ITS/eukaryome_ITS_dereplicated.classification
#
# Note: This script must be run from the project root directory.

suppressPackageStartupMessages({
  library(optparse)
  library(Biostrings)
  library(data.table)
  library(dplyr)
})

source("R/utils.R")

# ── Arguments ─────────────────────────────────────────────────────────────────

option_list <- list(
  make_option("--fasta_in",
              type = "character", metavar = "FILE",
              help = "Input FASTA file [required]"),
  make_option("--fasta_out",
              type = "character", metavar = "FILE",
              help = "Output FASTA file of dereplicated sequences [required]"),
  make_option("--classification_in",
              type = "character", metavar = "FILE",
              help = "Input tab-delimited classification file [required]"),
  make_option("--classification_out",
              type = "character", metavar = "FILE",
              help = "Output classification file with LCA taxonomy [required]")
)

opt <- parse_args(
  OptionParser(
    option_list = option_list,
    usage       = "%prog --fasta_in in.fa --fasta_out out.fa --classification_in in.tsv --classification_out out.tsv",
    description = "Dereplicate sequences and resolve taxonomy to the LCA."
  )
)

if (is.null(opt$fasta_in))           stop("--fasta_in is required.")
if (is.null(opt$fasta_out))          stop("--fasta_out is required.")
if (is.null(opt$classification_in))  stop("--classification_in is required.")
if (is.null(opt$classification_out)) stop("--classification_out is required.")
if (!file.exists(opt$fasta_in))          stop("File not found: ", opt$fasta_in)
if (!file.exists(opt$classification_in)) stop("File not found: ", opt$classification_in)

# ── Load inputs ───────────────────────────────────────────────────────────────

cat("Loading sequences from:", opt$fasta_in, "\n")
seqs <- readDNAStringSet(opt$fasta_in)
cat("Sequences loaded:", length(seqs), "\n")

cat("Loading classification from:", opt$classification_in, "\n")
cls <- fread(opt$classification_in)

# ── Dereplicate ───────────────────────────────────────────────────────────────

result <- dereplicate_lca(seqs, cls)

cat("\nDereplication summary\n")
cat("---------------------\n")
cat("Input sequences:      ", length(seqs), "\n")
cat("Unique sequences:     ", length(result$seqs), "\n")
cat("Removed duplicates:   ", length(seqs) - length(result$seqs), "\n")

# ── Write outputs ─────────────────────────────────────────────────────────────

writeXStringSet(result$seqs, opt$fasta_out)
cat("Dereplicated FASTA written to:", opt$fasta_out, "\n")

fwrite(result$classification, opt$classification_out, sep = "\t")
cat("LCA classification written to:", opt$classification_out, "\n")
