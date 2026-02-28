#!/usr/bin/env Rscript
# check_annotations.R — Detect and standardise infraspecific annotations
#
# Parses the species column of the classification file to extract infraspecific
# ranks (subsp., var., f., f.sp., subf.) and correct known formatting issues.
# The output file extends the classification with an 'infraspecific' column.
#
# Usage:
#   Rscript check_annotations.R \
#     --classification_in  data/full_ITS/eukaryome_ITS.classification \
#     --classification_out data/full_ITS/eukaryome_ITS.classification
#
# Note: This script must be run from the project root directory.

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(stringr)
})

# ── Arguments ─────────────────────────────────────────────────────────────────

option_list <- list(
  make_option("--classification_in",
              type = "character", metavar = "FILE",
              help = "Input tab-delimited classification file [required]"),
  make_option("--classification_out",
              type = "character", metavar = "FILE",
              help = "Output classification file with infraspecific column [required]")
)

opt <- parse_args(
  OptionParser(
    option_list = option_list,
    usage       = "%prog --classification_in in.tsv --classification_out out.tsv",
    description = paste(
      "Detect and standardise infraspecific annotations (subsp., var., f., etc.)",
      "in the species column of the classification file."
    )
  )
)

if (is.null(opt$classification_in))  stop("--classification_in is required.")
if (is.null(opt$classification_out)) stop("--classification_out is required.")
if (!file.exists(opt$classification_in)) stop("File not found: ", opt$classification_in)

class_in  <- opt$classification_in
class_out <- opt$classification_out

# ── Load classification ───────────────────────────────────────────────────────

cat("Loading classification from:", class_in, "\n")
taxa <- fread(class_in)

# ── Process species column ────────────────────────────────────────────────────

species <- taxa %>%
  select(id, species, genus) %>%
  mutate(
    # Fix known formatting issues
    species = str_replace_all(species, "\\+.1/2.", "e"),
    species = str_replace_all(species, "thaliana\\.x\\.arenosa\\.x\\.arenosa",
                                        "thaliana.x.arenosa"),
    species = str_replace_all(species, "metallica\\.x\\.Begonia\\.sanguinea",
                                        "metallica.x.sanguinea"),

    # Extract infraspecific annotation
    infraspecific = case_when(
      str_detect(species, "\\.subf\\.")  ~ str_extract(species, "\\.subf\\.[^\\s]+"),
      str_detect(species, "\\.f\\.sp\\.") ~ str_extract(species, "\\.f\\.sp\\.[^\\s]+"),
      str_detect(species, "\\.f\\.")     ~ str_extract(species, "\\.f\\.[^\\s]+"),
      str_detect(species, "\\.var\\.")   ~ str_extract(species, "\\.var\\.[^\\s]+"),
      str_detect(species, "\\.subsp\\.") ~ str_extract(species, "\\.subsp\\.[^\\s]+"),
      TRUE ~ NA_character_
    ),

    # Remove leading period from annotation
    infraspecific = str_remove(infraspecific, "^\\."),

    # Strip infraspecific suffix from species name
    species = case_when(
      !is.na(infraspecific) ~ str_remove(
        species,
        "\\.subf\\.[^\\s]+|\\.f\\.sp\\.[^\\s]+|\\.f\\.[^\\s]+|\\.var\\.[^\\s]+|\\.subsp\\.[^\\s]+"
      ),
      TRUE ~ species
    ),
    species       = str_trim(species),
    infraspecific = str_trim(infraspecific)
  )

# ── Summary ───────────────────────────────────────────────────────────────────

cat("\nInfraspecific annotations found:\n")
print(table(species$infraspecific, useNA = "ifany"))

# ── Detect remaining non-standard characters ──────────────────────────────────

unique_species <- unique(species$species)
species_split  <- tibble(species_name = unique_species) %>%
  filter(species_name != "unidentified") %>%
  mutate(
    genus_part   = str_extract(species_name, "^\\S+"),
    epithet_part = str_remove(species_name, "^\\S+\\s*")
  )

infraspecific_candidates <- species_split %>%
  filter(str_detect(epithet_part, "[^A-Za-z0-9_]")) %>%
  pull(species_name)

cat("\nSpecies with potential infraspecific annotations:", length(infraspecific_candidates), "\n")
if (length(infraspecific_candidates) > 0) {
  cat("First 20:\n")
  print(head(infraspecific_candidates, 20))
}

unique_abbr <- species_split %>%
  filter(species_name %in% infraspecific_candidates) %>%
  pull(epithet_part) %>%
  str_extract("\\.[^.]+(\\.([^.]+))?\\.") %>%
  .[!is.na(.)] %>% unique() %>% sort()
cat("\nUnique abbreviation patterns:\n")
print(unique_abbr)

# ── Manual corrections ────────────────────────────────────────────────────────

species <- species %>%
  mutate(
    species = case_when(
      species == "Oryza sativa.Indica"            ~ "Oryza sativa",
      species == "Oryza sativa.Japonica"           ~ "Oryza sativa",
      species == "Marchandiomyces nt59.1784"       ~ "Marchandiomyces unidentified",
      TRUE ~ species
    )
  )

# ── Rebuild classification ────────────────────────────────────────────────────

cat("\nReconstructing classification file...\n")

taxa_updated <- taxa %>%
  select(-species) %>%
  left_join(species %>% select(id, species, infraspecific), by = "id") %>%
  select(id, kingdom, phylum, class, order, family, genus, species, infraspecific)

fwrite(taxa_updated, class_out, sep = "\t")

cat("Written to:", class_out, "\n")
cat("Rows:", nrow(taxa_updated), "\n")
