#!/usr/bin/env Rscript
# test_load_classification_quoting.R — Regression test for the load_classification()
# quoting fix in R/predict.R.
#
# data.table::fwrite() (used by R/dereplicate_lca.R) writes a genuine
# empty-string value as the literal two characters "" -- to disambiguate it
# from NA, which it writes as nothing at all. Reading that back with
# quote = "" (as predict.R's load_classification() used to) disables quote
# parsing, so the literal "" is read back as two-character text rather than
# an empty string, and slips past is_identified() as if it were a real taxon
# name.
#
# This test builds a tiny classification file the same way dereplicate_lca.R
# would (an NA and a genuine "" value in the species column, written via
# fwrite), then confirms load_classification() (copied here, unchanged, from
# R/predict.R -- sourcing the whole file would trigger its CLI arg parsing)
# reads both back as blank / unidentified, not as literal text.
#
# Usage: Rscript tests/test_load_classification_quoting.R
# Note:  Must be run from the project root directory.

suppressPackageStartupMessages(library(data.table))

source("R/utils.R")

# ── load_classification(), copied verbatim from R/predict.R post-fix ─────────

load_classification <- function(class_file, ranks, higher_ranks, id_col = "id") {
  cls <- fread(class_file, sep = "\t", header = TRUE,
               fill = TRUE, na.strings = "", check.names = FALSE,
               data.table = FALSE)
  required <- unique(c(id_col, ranks, higher_ranks))
  missing  <- setdiff(required, colnames(cls))
  if (length(missing) > 0) {
    stop("Missing columns in classification file: ", paste(missing, collapse = ", "))
  }
  cls
}

# ── Build fixture the way dereplicate_lca.R's fwrite() actually would ────────

fixture_dir <- "tests/consolidate_cutoffs/tmp"
dir.create(fixture_dir, showWarnings = FALSE, recursive = TRUE)
fixture_path <- file.path(fixture_dir, "quoting_fixture.classification")

fixture <- data.table(
  id      = c("seq1", "seq2", "seq3"),
  kingdom = c("Fungi", "Fungi", "Fungi"),
  phylum  = c("Ascomycota", "Ascomycota", "Ascomycota"),
  class   = c("Sordariomycetes", "Sordariomycetes", "Sordariomycetes"),
  order   = c("Hypocreales", "Hypocreales", "Hypocreales"),
  family  = c("Clavicipitaceae", "Clavicipitaceae", "Clavicipitaceae"),
  genus   = c("Metarhizium", "Metarhizium", "Metarhizium"),
  species = c("Metarhizium anisopliae", NA_character_, "")
  # seq2: genuine NA        -> fwrite() writes nothing between the tabs
  # seq3: genuine ""        -> fwrite() writes the literal two chars ""
)

fwrite(fixture, fixture_path, sep = "\t")

# Sanity-check the fixture actually reproduces the bug's precondition: a
# literal "" in the file, not blank.
raw_lines <- readLines(fixture_path)
seq3_line <- raw_lines[grepl("^seq3", raw_lines)]
has_literal_quotes <- grepl('""', seq3_line, fixed = TRUE)

cat("Fixture file contents:\n")
cat(paste(raw_lines, collapse = "\n"), "\n\n")

n_pass <- 0L
n_fail <- 0L
check <- function(label, ok) {
  if (ok) { cat("PASS  ", label, "\n"); n_pass <<- n_pass + 1L }
  else    { cat("FAIL  ", label, "\n"); n_fail <<- n_fail + 1L }
}

check("fixture actually contains a literal \"\" for seq3 (precondition)", has_literal_quotes)

# ── The actual regression check ──────────────────────────────────────────────

cls <- load_classification(fixture_path, ranks = "species", higher_ranks = character(0))

check("seq1 species is identified (real name)", is_identified(cls$species[cls$id == "seq1"]))
check("seq2 (genuine NA) is NOT identified", !is_identified(cls$species[cls$id == "seq2"]))
check("seq3 (fwrite-escaped \"\") is NOT identified -- the actual regression",
      !is_identified(cls$species[cls$id == "seq3"]))
check("seq3 species reads back as a true empty string, not literal '\"\"'",
      identical(cls$species[cls$id == "seq3"], ""))

cat(sprintf("\nRESULT: %d passed, %d failed\n", n_pass, n_fail))
if (n_fail > 0) quit(status = 1)
