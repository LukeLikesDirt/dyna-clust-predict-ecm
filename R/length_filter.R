#!/usr/bin/env Rscript
# length_filter.R — Filter FASTA sequences by length and ambiguous nucleotides
#
# Retains sequences within [min_length, max_length]. By default also removes
# sequences containing IUPAC ambiguity codes (N R Y S W K M B D H V).
# Writes matching rows of the classification file.
#
# Usage:
#   Rscript length_filter.R \
#     --fasta_in              data/full_ITS/eukaryome_ITS.fasta \
#     --fasta_out             data/full_ITS/eukaryome_ITS_filtered.fasta \
#     --classification_in     data/full_ITS/eukaryome_ITS.classification \
#     --classification_out    data/full_ITS/eukaryome_ITS_filtered.classification \
#     --min_length 250 \
#     --max_length 1500 \
#     --exclude_ambiguous     # pass --no-exclude_ambiguous to skip
#
# Note: This script must be run from the project root directory.

suppressPackageStartupMessages({
  library(optparse)
  library(Biostrings)
  library(data.table)
})

source("R/utils.R")

# ── Arguments ─────────────────────────────────────────────────────────────────

option_list <- list(
  make_option("--fasta_in",
              type = "character", metavar = "FILE",
              help = "Input FASTA file [required]"),
  make_option("--fasta_out",
              type = "character", metavar = "FILE",
              help = "Output FASTA file of length-filtered sequences [required]"),
  make_option("--classification_in",
              type = "character", metavar = "FILE",
              help = "Input tab-delimited classification file [required]"),
  make_option("--classification_out",
              type = "character", metavar = "FILE",
              help = "Output classification file subset to surviving IDs [required]"),
  make_option("--min_length",
              type = "integer", default = 250L, metavar = "INT",
              help = "Minimum sequence length to retain [default: %default]"),
  make_option("--max_length",
              type = "integer", default = 1500L, metavar = "INT",
              help = "Maximum sequence length to retain [default: %default]"),
  make_option("--exclude_ambiguous",
              action = "store_true", default = TRUE,
              help = "Remove sequences with IUPAC ambiguity codes [default: %default]"),
  make_option("--no-exclude_ambiguous",
              action = "store_false", dest = "exclude_ambiguous",
              help = "Keep sequences with ambiguous nucleotides")
)

opt <- parse_args(
  OptionParser(
    option_list = option_list,
    usage       = "%prog --fasta_in in.fa --fasta_out out.fa --classification_in in.tsv --classification_out out.tsv",
    description = "Filter FASTA sequences by length and subset the classification file."
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

# ── Filter by length (and optionally ambiguous nucleotides) ───────────────────

n_after_length <- sum(Biostrings::width(seqs) >= opt$min_length &
                      Biostrings::width(seqs) <= opt$max_length)

seqs_filtered <- length_filter(seqs,
                               min_length         = opt$min_length,
                               max_length         = opt$max_length,
                               exclude_ambiguous  = opt$exclude_ambiguous)

cat("After length filter [", opt$min_length, "-", opt$max_length, "bp]:",
    n_after_length, "sequences\n")
if (opt$exclude_ambiguous) {
  cat("After ambiguous nucleotide filter:", length(seqs_filtered), "sequences\n")
}
cat("Total removed:", length(seqs) - length(seqs_filtered), "\n")

# ── Subset classification to surviving IDs ────────────────────────────────────

cls_filtered <- cls[cls$id %in% names(seqs_filtered), ]
cat("Classification rows after filter:", nrow(cls_filtered), "\n")

# ── Write outputs ─────────────────────────────────────────────────────────────

writeXStringSet(seqs_filtered, opt$fasta_out)
cat("Filtered FASTA written to:", opt$fasta_out, "\n")

fwrite(cls_filtered, opt$classification_out, sep = "\t")
cat("Filtered classification written to:", opt$classification_out, "\n")
