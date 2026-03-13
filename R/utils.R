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

# ── Length filtering ──────────────────────────────────────────────────────────
#
# length_filter(seqs, min_length = 250, max_length = 1500, exclude_ambiguous = TRUE)
#
# Filters a DNAStringSet to retain only sequences within [min_length, max_length].
# When exclude_ambiguous = TRUE (default), sequences containing any IUPAC
# ambiguity codes (N R Y S W K M B D H V) are also removed.
# Returns the filtered DNAStringSet.

length_filter <- function(seqs, min_length = 250, max_length = 1500,
                          exclude_ambiguous = TRUE) {
  lens <- Biostrings::width(seqs)
  keep <- lens >= min_length & lens <= max_length
  seqs <- seqs[keep]

  if (exclude_ambiguous) {
    ambig_counts <- rowSums(
      Biostrings::letterFrequency(seqs, letters = c("N","R","Y","S","W","K","M","B","D","H","V"))
    )
    seqs <- seqs[ambig_counts == 0]
  }

  seqs
}

# ── Dereplication with LCA taxonomy ──────────────────────────────────────────
#
# dereplicate_lca(seqs, cls)
#
# Collapses identical sequences and resolves taxonomy to the lowest common
# ancestor (LCA). Sequences with identical strings are grouped; within each
# group the first ID is kept as representative and ranks are resolved top-down,
# stopping at the first rank where two or more distinct non-"unidentified"
# values exist.
#
# Arguments:
#   seqs  DNAStringSet with named sequences.
#   cls   data.frame / data.table with columns: id, kingdom, phylum, class,
#         order, family, genus, species.
#
# Returns a named list:
#   $seqs            DNAStringSet of representative sequences.
#   $classification  data.frame with id + taxonomy columns.

dereplicate_lca <- function(seqs, cls) {

  taxonomy_ranks <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")

  seq_strings <- as.character(seqs)
  derep_table <- data.frame(
    id       = names(seqs),
    sequence = seq_strings,
    stringsAsFactors = FALSE
  )

  derep_tax <- derep_table %>%
    dplyr::left_join(as.data.frame(cls), by = "id")

  resolve_lca <- function(df) {
    resolved <- list()
    for (rank in taxonomy_ranks) {
      vals <- unique(df[[rank]])
      vals <- vals[!vals %in% c("unidentified", NA)]
      if (length(vals) == 1) {
        resolved[[rank]] <- vals
      } else {
        resolved[[rank]] <- NA
        break
      }
    }
    as.data.frame(resolved)
  }

  resolved_tax <- derep_tax %>%
    dplyr::group_by(sequence) %>%
    dplyr::group_modify(~ resolve_lca(.x)) %>%
    dplyr::ungroup()

  rep_ids <- derep_table %>%
    dplyr::group_by(sequence) %>%
    dplyr::summarise(id = dplyr::first(id), .groups = "drop")

  resolved_tax <- resolved_tax %>%
    dplyr::left_join(rep_ids, by = "sequence")

  rep_seqs        <- seqs[rep_ids$id]
  names(rep_seqs) <- rep_ids$id

  classification_out <- resolved_tax %>%
    dplyr::select(id, dplyr::all_of(taxonomy_ranks)) %>%
    dplyr::arrange(id)

  list(seqs = rep_seqs, classification = as.data.frame(classification_out))
}
