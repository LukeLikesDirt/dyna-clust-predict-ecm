#!/usr/bin/env Rscript
# reformat.R — Reformat EUKARYOME FASTA headers and extract taxonomy
#
# Usage:
#   Rscript reformat.R --fasta_in raw.fasta --fasta_out clean.fasta \
#                      --classification_out taxonomy.tsv
#
# Note: This script must be run from the project root directory.

required_packages <- c("optparse", "dplyr", "readr", "stringr", "tidyr", "furrr", "future")
missing_packages  <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing packages: ", paste(missing_packages, collapse = ", "),
       "\nInstall with: install.packages(c('", paste(missing_packages, collapse = "', '"), "'))")
}

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(furrr)
  library(future)
})

options(future.globals.maxSize = 2000 * 1024^2)

# ── Arguments ─────────────────────────────────────────────────────────────────

option_list <- list(
  make_option("--fasta_in",
              type = "character", metavar = "FILE",
              help = "Input raw FASTA file [required]"),
  make_option("--fasta_out",
              type = "character", metavar = "FILE",
              help = "Output reformatted FASTA file [required]"),
  make_option("--classification_out",
              type = "character", metavar = "FILE",
              help = "Output tab-delimited classification file [required]")
)

opt <- parse_args(
  OptionParser(
    option_list = option_list,
    usage       = "%prog --fasta_in raw.fasta --fasta_out clean.fasta --classification_out taxonomy.tsv",
    description = paste(
      "Reformat EUKARYOME FASTA headers for dyna-clust-predict.",
      "Extracts taxonomic ranks, standardises names, deduplicates by ID."
    )
  )
)

if (is.null(opt$fasta_in))           stop("--fasta_in is required.")
if (is.null(opt$fasta_out))          stop("--fasta_out is required.")
if (is.null(opt$classification_out)) stop("--classification_out is required.")
if (!file.exists(opt$fasta_in))      stop("Input FASTA not found: ", opt$fasta_in)

fasta_in           <- opt$fasta_in
fasta_out          <- opt$fasta_out
classification_out <- opt$classification_out

for (d in unique(dirname(c(fasta_out, classification_out)))) {
  if (!dir.exists(d) && d != ".") dir.create(d, recursive = TRUE)
}

# ── Parallel FASTA reading ────────────────────────────────────────────────────

cat("Reading FASTA file:", fasta_in, "\n")
fasta_lines <- readLines(fasta_in)

header_indices <- which(grepl("^>", fasta_lines))
cat("Processing", length(header_indices), "sequences in parallel...\n")

plan(multicore, workers = availableCores() - 1)

fasta_df <- future_map_dfr(seq_along(header_indices), function(i) {
  h   <- header_indices[i]
  end <- if (i < length(header_indices)) header_indices[i + 1] - 1L else length(fasta_lines)
  tibble(
    header   = fasta_lines[h],
    sequence = paste(fasta_lines[(h + 1):end], collapse = "")
  )
}, .progress = TRUE, .options = furrr_options(seed = TRUE))

plan(sequential)

# ── Taxonomy processing ───────────────────────────────────────────────────────

cat("\nProcessing taxonomic classifications...\n")

classification_df <- fasta_df %>%
  mutate(
    id             = str_extract(header, "(?<=^>)[^;]+"),
    kingdom        = str_remove(str_extract(header, "k__([^;]+)"), "k__"),
    phylum         = str_remove(str_extract(header, "p__([^;]+)"), "p__"),
    class          = str_remove(str_extract(header, "c__([^;]+)"), "c__"),
    order          = str_remove(str_extract(header, "o__([^;]+)"), "o__"),
    family         = str_remove(str_extract(header, "f__([^;]+)"), "f__"),
    genus          = str_remove(str_extract(header, "g__([^;]+)"), "g__"),
    species_epithet = str_remove(str_extract(header, "s__([^;]+)$"), "s__")
  ) %>%
  # Standardise delimiters: replace "-" with "."
  mutate(across(c(kingdom, phylum, class, order, family, genus, species_epithet),
                ~ str_replace_all(.x, "-", "."))) %>%
  # Replace "unclassified" with "unidentified"
  mutate(across(c(kingdom, phylum, class, order, family, genus, species_epithet),
                ~ ifelse(.x == "unclassified", "unidentified", .x))) %>%
  # Remove tentative IDs (cf., nr., aff., nov.inval.) and EUKARYOME predictions
  mutate(
    across(
      c(kingdom, phylum, class, order, family, genus, species_epithet),
      ~ ifelse(str_detect(.x, "cf\\.|nr\\.|aff\\.|nov\\.inval\\.|\\.reg"),
               "unidentified", .x)
    ),
    phylum = ifelse(grepl("\\.phy", phylum), "unidentified", phylum),
    class  = ifelse(grepl("\\.cl",  class),  "unidentified", class),
    order  = ifelse(grepl("\\.ord", order),  "unidentified", order),
    family = ifelse(grepl("\\.fam", family), "unidentified", family),
    genus  = ifelse(grepl("\\.gen", genus),  "unidentified", genus)
  ) %>%
  # Parenthetical names: replace ( ) with __ for downstream parsing
  mutate(
    across(
      c(phylum, class, order, family, genus),
      ~ ifelse(grepl("\\(", .x),
               str_replace_all(str_replace_all(.x, "\\(", "__"), "\\)", "__"),
               .x)
    )
  ) %>%
  # Construct species names
  mutate(
    species = case_when(
      is.na(species_epithet) | species_epithet == "" | species_epithet == "unidentified" ~ "unidentified",
      is.na(genus)           | genus == ""           | genus == "unidentified"           ~ "unidentified",
      TRUE ~ paste(genus, species_epithet, sep = " ")
    )
  ) %>%
  # Taxonomic completeness score for deduplication
  mutate(
    tax_score =
      (!is.na(kingdom) & kingdom != "" & kingdom != "unidentified") +
      (!is.na(phylum)  & phylum  != "" & phylum  != "unidentified") +
      (!is.na(class)   & class   != "" & class   != "unidentified") +
      (!is.na(order)   & order   != "" & order   != "unidentified") +
      (!is.na(family)  & family  != "" & family  != "unidentified") +
      (!is.na(genus)   & genus   != "" & genus   != "unidentified") +
      (!is.na(species) & species != "" & species != "unidentified")
  ) %>%
  group_by(id) %>%
  arrange(desc(tax_score), id) %>%
  dplyr::slice(1) %>%
  ungroup() %>%
  mutate(across(everything(), ~ ifelse(is.na(.x), "", .x)))

if (any(duplicated(classification_df$id))) {
  dup_ids <- classification_df$id[duplicated(classification_df$id)]
  stop("Duplicate IDs found: ", paste(head(dup_ids, 10), collapse = ", "))
}

# ── Write outputs ─────────────────────────────────────────────────────────────

cat("Writing output files...\n")

fasta_content <- classification_df %>%
  mutate(fasta_header = paste0(">", id)) %>%
  select(fasta_header, sequence) %>%
  pivot_longer(everything(), names_to = NULL, values_to = "line") %>%
  pull(line)
writeLines(fasta_content, fasta_out)

classification_df %>%
  select(id, kingdom, phylum, class, order, family, genus, species) %>%
  write_tsv(classification_out)

# ── Summary ───────────────────────────────────────────────────────────────────

cat("\nDone.\n")
cat("  Input:           ", fasta_in, "\n")
cat("  Output FASTA:    ", fasta_out, "\n")
cat("  Classification:  ", classification_out, "\n")
cat("  Sequences:       ", nrow(classification_df), "\n\n")
for (rank in c("phylum", "class", "order", "family", "genus", "species")) {
  n_id   <- sum(classification_df[[rank]] != "unidentified")
  n_unid <- sum(classification_df[[rank]] == "unidentified")
  cat(sprintf("  %-10s identified: %d  unidentified: %d\n", rank, n_id, n_unid))
}
