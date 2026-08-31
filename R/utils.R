# utils.R — Shared utility functions for dyna-clust-predict
#
# Source this file at the top of any script that needs taxonomy filtering:
#   source("R/utils.R")

# ── Rank hierarchy metadata ───────────────────────────────────────────────────
#
# Shared by R/subset.R and R/consolidate_cutoffs.R so the taxonomic hierarchy
# is defined exactly once.
#
# rank_hierarchy   Coarsest -> finest.
# rank_abbr        Filename abbreviation used for <target>_pred_id_<abbr>.txt.
# parent_ranks_map Target rank -> valid parent ranks, coarsest -> finest is
#                  reversed here (finest parent first) since nested_prediction_filter()
#                  and consolidate_cutoffs.R both walk parents fine-to-coarse.

rank_hierarchy <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")

rank_abbr <- c(
  kingdom = "kng", phylum = "phy", class = "cls", order = "ord",
  family  = "fam", genus  = "gen", species = "spe"
)

parent_ranks_map <- list(
  species = c("genus", "family", "order", "class", "phylum", "kingdom"),
  genus   = c("family", "order", "class", "phylum", "kingdom"),
  family  = c("order", "class", "phylum", "kingdom"),
  order   = c("class", "phylum", "kingdom"),
  class   = c("phylum", "kingdom"),
  phylum  = c("kingdom"),
  kingdom = character(0)
)

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
#   - starts with "_" (EUKARYOME internal/artifact markers, e.g.
#     _mitochondrion, _nucleomorph, _plastid, _Archaea, _Bacteria,
#     _pseudogene, _pseudogene.Fungi -- no legitimate taxon name at any rank
#     starts with an underscore, so this single rule covers the whole family
#     of markers; the explicit alternation below is kept for documentation)

is_identified <- function(x) {
  !is.na(x) &
  nzchar(x) &
  !grepl("^unidentified$|^unclassified$|^uncultured", x, ignore.case = TRUE) &
  !grepl("incertae[ _.-]sedis",                       x, ignore.case = TRUE) &
  !grepl("[ _.-]sp\\.",                               x, ignore.case = TRUE) &
  !grepl("Unispike1|Unispike2|Unispike3",             x, ignore.case = TRUE) &
  !grepl("Archaea|Bacteria",                          x, ignore.case = TRUE) &
  !grepl("mitochondrion|nucleomorph|plastid",         x, ignore.case = TRUE) &
  !grepl("^_",                                        x)
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
      vals <- vals[!vals %in% c("unidentified", "unclassified", NA, "")]
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
    dplyr::arrange(id) %>%
    dplyr::mutate(dplyr::across(dplyr::all_of(taxonomy_ranks),
                                ~ dplyr::if_else(is.na(.x) | .x == "", "unclassified", .x)))

  list(seqs = rep_seqs, classification = as.data.frame(classification_out))
}
