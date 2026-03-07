# utils.R — Shared utility functions for dyna-clust-predict
#
# Source this file at the top of any script that needs taxonomy filtering:
#   source("R/utils.R")

# ── Taxonomy filtering ────────────────────────────────────────────────────────
#
# is_identified(x)
#
# Returns a logical vector: TRUE where x represents a valid, resolved
# taxonomic name; FALSE for missing, empty, or ambiguous entries.
#
# Patterns treated as unidentified (case-insensitive):
#   - NA or empty string
#   - "unidentified" or "unclassified" (exact, full word)
#   - starts with "uncultured"
#   - contains "incertae sedis" (with space, underscore, dot, or dash separator)
#   - ends with " sp.", "_sp.", ".sp.", or "-sp." (species placeholder)

is_identified <- function(x) {
  !is.na(x) &
  nzchar(x) &
  !grepl("^unidentified$|^unclassified$|^uncultured", x, ignore.case = TRUE) &
  !grepl("incertae[ _.-]sedis",                       x, ignore.case = TRUE) &
  !grepl("[ _.-]sp\\.",                               x, ignore.case = TRUE) &
  !grepl("Unispike1|Unispike2|Unispike3",             x, ignore.case = TRUE) &
  !grepl("Archaea|Bacteria",                          x, ignore.case = TRUE) &
  !grepl("mitochondrion|nucleomorph|plastid",         x, ignore.case = TRUE) 
}
