#!/usr/bin/env Rscript
# run_test.R — Correctness test for R/reformat.R and the is_identified()
# artifact-marker fix in R/utils.R.
#
# Two sections:
#   1. Direct is_identified() checks (fast, no subprocess) -- confirms the
#      new leading-underscore rule rejects EUKARYOME's artifact/pseudogene
#      kingdom markers without disturbing the pre-existing explicit checks.
#   2. A full R/reformat.R subprocess run against a small synthetic FASTA
#      fixture (exactly how the real pipeline invokes it), asserting the
#      resulting classification + homonym manifest cell-by-cell.
#
# Fixture design -- one scenario per synthetic accession:
#   s1        control: ordinary genus/species, no collision, no placeholder
#   s2 (Fusarium) shares the same genus as s1 deliberately, to also serve as
#             the "non-colliding genus left untouched" check
#   s3/s4     Alsophila under Metazoa and Viridiplantae -- a genuine
#             cross-kingdom homonym; both genus AND species must pick up the
#             "(Kingdom)" suffix
#   s5        Galeropsis(Fungi) -- an upstream-supplied suffix already
#             present in the raw header; must pass through unchanged (guards
#             against double-suffixing) and must NOT appear in the manifest
#   s6        (Aleuria) -- fully-wrapped genus, well-formed parens
#   s7        (Candida] -- fully-wrapped genus, MISMATCHED closing bracket
#   s8/s9     Entrophospora under _pseudogene and under real Fungi -- must
#             NOT be treated as a homonym once _pseudogene is correctly
#             rejected by is_identified()
#   s10       Densosporales.fam02 / Densosporales.fam02.gen01 -- Tedersoo
#             et al. (2024) placeholder codes, must survive with "." -> "_"
#
# Usage: Rscript tests/reformat/run_test.R
# Note:  Must be run from the project root directory.

suppressPackageStartupMessages(library(data.table))

source("R/utils.R")

n_pass <- 0L
n_fail <- 0L
check <- function(label, ok) {
  if (isTRUE(ok)) { cat("PASS  ", label, "\n"); n_pass <<- n_pass + 1L }
  else            { cat("FAIL  ", label, "\n"); n_fail <<- n_fail + 1L }
}

# ══════════════════════════════════════════════════════════════════════════════
# Section 1: is_identified() artifact-marker checks
# ══════════════════════════════════════════════════════════════════════════════

cat("══════════════════════════════════════════════════════════════\n")
cat("Section 1: is_identified() artifact-marker rejection\n")
cat("══════════════════════════════════════════════════════════════\n")

check("_pseudogene rejected (the actual regression)", !is_identified("_pseudogene"))
check("_pseudogene.Fungi rejected", !is_identified("_pseudogene.Fungi"))
check("_mitochondrion still rejected (pre-existing)", !is_identified("_mitochondrion"))
check("_nucleomorph still rejected (pre-existing)", !is_identified("_nucleomorph"))
check("_plastid still rejected (pre-existing)", !is_identified("_plastid"))
check("_Archaea still rejected (pre-existing)", !is_identified("_Archaea"))
check("_Bacteria still rejected (pre-existing)", !is_identified("_Bacteria"))
check("real kingdom Fungi still accepted", is_identified("Fungi"))
check("real kingdom Viridiplantae still accepted", is_identified("Viridiplantae"))
check("real genus with underscore placeholder still accepted (not a leading underscore)",
      is_identified("Densosporales_fam02"))

# ══════════════════════════════════════════════════════════════════════════════
# Section 2: full R/reformat.R subprocess run
# ══════════════════════════════════════════════════════════════════════════════

cat("\n══════════════════════════════════════════════════════════════\n")
cat("Section 2: R/reformat.R end-to-end\n")
cat("══════════════════════════════════════════════════════════════\n")

fixture_dir <- "tests/reformat/tmp"
dir.create(fixture_dir, showWarnings = FALSE, recursive = TRUE)

fasta_in     <- file.path(fixture_dir, "raw.fasta")
fasta_out    <- file.path(fixture_dir, "reformatted.fasta")
class_out    <- file.path(fixture_dir, "reformatted.classification")
manifest_out <- file.path(fixture_dir, "homonym_manifest.txt")

header <- function(id, kingdom, phylum, class, order, family, genus, species) {
  sprintf(">%s;k__%s;p__%s;c__%s;o__%s;f__%s;g__%s;s__%s",
          id, kingdom, phylum, class, order, family, genus, species)
}

fasta_lines <- c(
  header("s1", "Fungi", "Ascomycota", "Sordariomycetes", "Hypocreales",
         "Nectriaceae", "Fusarium", "oxysporum"),
  "ACGTACGTACGT",

  header("s2", "Fungi", "Ascomycota", "Sordariomycetes", "Hypocreales",
         "Nectriaceae", "Fusarium", "solani"),
  "ACGTACGTACGT",

  header("s3", "Metazoa", "Arthropoda", "Insecta", "Lepidoptera",
         "Geometridae", "Alsophila", "pometaria"),
  "ACGTACGTACGT",

  header("s4", "Viridiplantae", "Tracheophyta", "Polypodiopsida", "Cyatheales",
         "Cyatheaceae", "Alsophila", "spinulosa"),
  "ACGTACGTACGT",

  header("s5", "Fungi", "Basidiomycota", "Agaricomycetes", "Agaricales",
         "Bolbitiaceae", "Galeropsis(Fungi)", "capensis"),
  "ACGTACGTACGT",

  header("s6", "Fungi", "Ascomycota", "Pezizomycetes", "Pezizales",
         "Pyronemataceae", "(Aleuria)", "aurantia"),
  "ACGTACGTACGT",

  header("s7", "Fungi", "Ascomycota", "Saccharomycetes", "Serinales",
         "Debaryomycetaceae", "(Candida]", "glabrata"),
  "ACGTACGTACGT",

  header("s8", "_pseudogene", "Glomeromycota", "Glomeromycetes",
         "Entrophosporales", "Entrophosporaceae", "Entrophospora", "sp01"),
  "ACGTACGTACGT",

  header("s9", "Fungi", "Glomeromycota", "Glomeromycetes",
         "Entrophosporales", "Entrophosporaceae", "Entrophospora", "infrequens"),
  "ACGTACGTACGT",

  header("s10", "Fungi", "Glomeromycota", "Glomeromycetes", "Densosporales",
         "Densosporales.fam02", "Densosporales.fam02.gen01", "indet"),
  "ACGTACGTACGT"
)

writeLines(fasta_lines, fasta_in)

rc <- system2(
  "Rscript",
  args = c("R/reformat.R",
           "--fasta_in", fasta_in,
           "--fasta_out", fasta_out,
           "--classification_out", class_out,
           "--manifest_out", manifest_out),
  stdout = TRUE, stderr = TRUE
)
cat(paste(rc, collapse = "\n"), "\n\n")

exit_status <- attr(rc, "status")
if (!is.null(exit_status) && exit_status != 0) {
  cat("*** FAIL: R/reformat.R exited with status", exit_status, "***\n")
  quit(status = 1)
}
if (!file.exists(class_out)) {
  cat("*** FAIL: classification output was not created:", class_out, "***\n")
  quit(status = 1)
}

cls <- fread(class_out, sep = "\t", header = TRUE)

get_row <- function(seq_id) {
  # Parameter named seq_id (not "id") to avoid colliding with the `id`
  # column inside data.table's [ NSE -- a bare `id` on the right of `==`
  # would resolve to the column itself, matching every row against itself.
  hit <- cls[id == seq_id]
  if (nrow(hit) == 0) NULL else hit[1]
}

check_field <- function(label, seq_id, field, expected) {
  row <- get_row(seq_id)
  if (is.null(row)) {
    cat(sprintf("FAIL  %-60s  row not found (%s)\n", label, seq_id))
    n_fail <<- n_fail + 1L
    return(invisible())
  }
  actual <- row[[field]]
  ok <- identical(actual, expected)
  if (ok) cat(sprintf("PASS  %-60s  %s = %s\n", label, field, actual))
  else    cat(sprintf("FAIL  %-60s  %s = %s (expected %s)\n", label, field, actual, expected))
  if (ok) n_pass <<- n_pass + 1L else n_fail <<- n_fail + 1L
}

cat("\n-- control: non-colliding genus untouched --\n")
check_field("s1.genus unchanged (Fusarium, shared but same kingdom -> no collision)",
            "s1", "genus", "Fusarium")
check_field("s2.genus unchanged (same as s1)", "s2", "genus", "Fusarium")
check_field("s1.species constructed normally", "s1", "species", "Fusarium oxysporum")

cat("\n-- genuine cross-kingdom homonym: disambiguated in genus AND species --\n")
check_field("s3.genus disambiguated (Metazoa)", "s3", "genus", "Alsophila(Metazoa)")
check_field("s3.species inherits disambiguated genus",
            "s3", "species", "Alsophila(Metazoa) pometaria")
check_field("s4.genus disambiguated (Viridiplantae)", "s4", "genus", "Alsophila(Viridiplantae)")
check_field("s4.species inherits disambiguated genus",
            "s4", "species", "Alsophila(Viridiplantae) spinulosa")

cat("\n-- upstream-supplied suffix: passed through, not double-wrapped --\n")
check_field("s5.genus unchanged (already suffixed upstream)",
            "s5", "genus", "Galeropsis(Fungi)")
check_field("s5.species unchanged prefix", "s5", "species", "Galeropsis(Fungi) capensis")

cat("\n-- fully-wrapped genus stripped (well-formed and mismatched brackets) --\n")
check_field("s6.genus (Aleuria) -> Aleuria", "s6", "genus", "Aleuria")
check_field("s7.genus (Candida] -> Candida (mismatched bracket)", "s7", "genus", "Candida")

cat("\n-- _pseudogene does not manufacture a false homonym --\n")
check_field("s8.genus unchanged (own kingdom _pseudogene is unidentified)",
            "s8", "genus", "Entrophospora")
check_field("s9.genus unchanged (only 1 identified kingdom: Fungi)",
            "s9", "genus", "Entrophospora")

cat("\n-- Tedersoo et al. (2024) placeholder codes preserved --\n")
check_field("s10.family placeholder preserved, '.' -> '_'",
            "s10", "family", "Densosporales_fam02")
check_field("s10.genus placeholder preserved, '.' -> '_'",
            "s10", "genus", "Densosporales_fam02_gen01")

cat("\n-- homonym manifest --\n")
if (!file.exists(manifest_out)) {
  cat("FAIL   manifest file was not created\n")
  n_fail <- n_fail + 1L
} else {
  manifest <- fread(manifest_out, sep = "\t", header = TRUE)
  check("manifest contains exactly 2 rows (Alsophila x 2 kingdoms)", nrow(manifest) == 2L)
  check("manifest does not mention Galeropsis (already disambiguated upstream)",
        !("Galeropsis(Fungi)" %in% manifest$original_genus) &&
        !("Galeropsis" %in% manifest$original_genus))
  check("manifest does not mention Entrophospora (_pseudogene excluded)",
        !("Entrophospora" %in% manifest$original_genus))
  check("manifest does not mention Fusarium (no collision)",
        !("Fusarium" %in% manifest$original_genus))
  also <- manifest[manifest$original_genus == "Alsophila"]
  check("manifest has an Alsophila/Metazoa row",
        any(also$kingdom == "Metazoa" & also$disambiguated_genus == "Alsophila(Metazoa)"))
  check("manifest has an Alsophila/Viridiplantae row",
        any(also$kingdom == "Viridiplantae" & also$disambiguated_genus == "Alsophila(Viridiplantae)"))
}

cat(sprintf("\nRESULT: %d passed, %d failed\n", n_pass, n_fail))
if (n_fail > 0) quit(status = 1)
