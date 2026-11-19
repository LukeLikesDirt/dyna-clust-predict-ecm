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

# ── Lineage canonicalisation ──────────────────────────────────────────────────
#
# EUKARYOME is not a tree: some names sit under more than one parent, e.g.
# Cladosporiaceae under Cladosporiales (564 ITS2 records) and Dothideales (344),
# Candida under 6 families, Triscelophorus under two spellings of its family,
# order Talbotiomycetales under class Talbotiomycetes and under a class named
# Talbotiomycetales. Left alone, consolidate_cutoffs.R warns ("inconsistent
# parent ... modal parent used for lineage lookup") and walks such a taxon's
# ancestors by its modal parent; dereplicate_lca.R truncates the LCA of an
# identical-sequence group at the conflicting rank; and a downstream classifier
# sees one name under two parents.
#
# canonicalise_lineage() gives every named label exactly one parent, so the
# classification file is a tree from the point it is generated (step 01).
# Same rule as canonicalise_lineage() in the ITS2 clustering pipeline that
# consumes these files (molly/bioinformatics/its2/utils/taxonomy.R), so that
# pipeline's own pass finds nothing left to change.
#
# is_missing_label(x)  NA, "", "unidentified" or "unclassified": no label at
#                      that rank, so nothing to canonicalise.
# is_incertae(x)       an incertae sedis label (a placement, but not a name).

is_missing_label <- function(x) {
  is.na(x) | !nzchar(x) | x %in% c("unidentified", "unclassified")
}

is_incertae <- function(x) !is.na(x) & grepl("incertae[ _.-]sedis", x, ignore.case = TRUE)

# find_nesting_violations(cls, ranks = rank_hierarchy)
#
# Labels that appear under more than one parent label (named child ranks only;
# a missing parent label counts as a parent). Returns a data.table
# (rank, label, n_parents); zero rows means the table is a tree.
# Requires data.table to be attached by the calling script.

find_nesting_violations <- function(cls, ranks = rank_hierarchy) {
  out <- vector("list", length(ranks) - 1L)
  for (k in 2:length(ranks)) {
    x <- unique(data.table(child = as.character(cls[[ranks[k]]]),
                           parent = as.character(cls[[ranks[k - 1L]]])))
    x <- x[!is_missing_label(child)]
    v <- x[, .(n_parents = .N), by = child][n_parents > 1L]
    out[[k - 1L]] <- data.table(rank = rep(ranks[k], nrow(v)), label = v$child, n_parents = v$n_parents)
  }
  rbindlist(out)
}

# canonicalise_lineage(cls, ranks = rank_hierarchy)
#
#   cls  data.frame / data.table with `id`, the rank columns and anything else
#        (kept as is, in the same row order).
#
# For each child rank, species up to phylum, the parent of a named label is its
# modal parent among identified parents (is_identified()), or among incertae
# sedis parents if it has no identified one, or among missing labels if it has
# neither. An identified parent therefore beats an incertae one however few
# records carry it, and a named label whose parent is missing in some records
# gets its identified parent in all of them. Ties go to the alphabetically
# first parent (C locale). The vote is taken on the ORIGINAL table at every
# rank, then applied bottom-up to every record carrying that label, so a
# record's whole lineage becomes the chain of chosen parents above its finest
# named label. Labels missing at a rank are left alone.
#
# Records under an EUKARYOME marker kingdom (a kingdom starting with "_":
# _mitochondrion, _pseudogene, _plastid ...) neither vote nor change, so a
# marker is never re-parented into a real kingdom (EUKARYOME files 42 ITS
# records under kingdom _pseudogene, phylum Glomeromycota).
#
# Returns list(classification, manifest, n_rows_changed, n_orphans):
#   classification  same rows, columns and order as `cls`, now a tree
#   manifest        one row per label that had more than one parent: rank,
#                   child, chosen_parent, n_parents, n_records,
#                   n_records_changed, parents ("A(786); B(4)")
#   n_rows_changed  records whose lineage differs from the input
#   n_orphans       per rank, records with a named label under a missing
#                   parent label that has no identified parent anywhere
#                   (reported, not changed)
# Stops if a label still has two parents afterwards, or if a record would move
# from one identified kingdom to another (a cross-kingdom homonym that
# reformat.R's "(Kingdom)" suffix does not cover).
# Requires data.table to be attached by the calling script.

canonicalise_lineage <- function(cls, ranks = rank_hierarchy) {
  cls <- copy(as.data.table(cls))
  if (!all(c("id", ranks) %in% names(cls))) {
    stop("cls needs columns: id, ", paste(ranks, collapse = ", "))
  }
  if (any(c("uid__", "marker__") %in% names(cls))) stop("cls may not have columns uid__ or marker__")
  for (r in ranks) set(cls, j = r, value = as.character(cls[[r]]))

  k1 <- cls[[ranks[1L]]]
  set(cls, j = "marker__", value = !is.na(k1) & startsWith(k1, "_"))
  cls[, uid__ := .GRP, by = c(ranks, "marker__")]
  u <- cls[, c(list(N = .N, marker = marker__[1L]), lapply(.SD, function(z) z[1L])),
           by = uid__, .SDcols = ranks]
  setorder(u, uid__)

  maps     <- setNames(vector("list", length(ranks)), ranks)
  manifest <- list()
  for (k in length(ranks):2L) {
    child_rank  <- ranks[k]
    parent_rank <- ranks[k - 1L]
    v <- data.table(child = u[[child_rank]], parent = u[[parent_rank]], N = u$N)[!u$marker]
    v <- v[!is_missing_label(child), .(N = sum(N)), by = .(child, parent)]
    if (!nrow(v)) {
      maps[[child_rank]] <- data.table(child = character(0), parent = character(0))
      next
    }
    v[, preference := fifelse(is_identified(parent), 1L, fifelse(is_incertae(parent), 2L, 3L))]
    v[, best_pref := min(preference), by = child]
    cand <- v[preference == best_pref]
    setorder(cand, child, -N, parent)
    chosen <- cand[, .(parent = parent[1L]), by = child]
    maps[[child_rank]] <- chosen

    tot   <- v[, .(n_parents = .N, n_records = sum(N)), by = child]
    multi <- tot[n_parents > 1L]
    if (nrow(multi)) {
      listed <- v[child %in% multi$child][order(child, -N, parent),
                   .(parents = paste0(parent, "(", N, ")", collapse = "; ")), by = child]
      n_kept <- merge(chosen, v, by = c("child", "parent"))[, .(child, n_kept = N)]
      m <- Reduce(function(a, b) merge(a, b, by = "child"), list(multi, chosen, listed, n_kept))
      m[, `:=`(rank = child_rank, chosen_parent = parent, n_records_changed = n_records - n_kept)]
      manifest[[child_rank]] <- m[, .(rank, child, chosen_parent, n_parents, n_records,
                                      n_records_changed, parents)]
    }
  }

  new <- copy(u)
  for (k in length(ranks):2L) {
    child_rank  <- ranks[k]
    parent_rank <- ranks[k - 1L]
    ch    <- new[[child_rank]]
    named <- which(!is_missing_label(ch) & !new$marker)
    if (!length(named)) next
    mp  <- maps[[child_rank]]
    idx <- match(ch[named], mp$child)
    if (anyNA(idx)) stop("internal error: a named ", child_rank, " label has no vote")
    set(new, i = named, j = parent_rank, value = mp$parent[idx])
  }

  moved <- is_identified(u[[ranks[1L]]]) & is_identified(new[[ranks[1L]]]) &
           u[[ranks[1L]]] != new[[ranks[1L]]]
  if (any(moved)) {
    i <- which(moved)[1L]
    stop(sum(u$N[moved]), " record(s) would move to another ", ranks[1L], ", e.g. from ",
         u[[ranks[1L]]][i], " to ", new[[ranks[1L]]][i], " (lineage ",
         paste(unlist(u[i, ranks, with = FALSE]), collapse = ";"),
         "): a cross-kingdom homonym that needs a \"(Kingdom)\" suffix in reformat.R")
  }

  changed <- Reduce(`|`, lapply(ranks, function(r) new[[r]] != u[[r]]))
  n_rows_changed <- sum(u$N[changed])

  viol <- find_nesting_violations(new[marker == FALSE], ranks)
  if (nrow(viol)) {
    stop(nrow(viol), " label(s) still have more than one parent after canonicalising (first: ",
         viol$rank[1L], " / ", viol$label[1L], ")")
  }
  n_orphans <- vapply(2:length(ranks), function(k) {
    keep <- !new$marker & !is_missing_label(new[[ranks[k]]]) & is_missing_label(new[[ranks[k - 1L]]])
    as.integer(sum(new$N[keep]))
  }, integer(1))
  names(n_orphans) <- ranks[-1L]

  out <- cls[, setdiff(names(cls), c("uid__", "marker__")), with = FALSE]
  for (r in ranks) set(out, j = r, value = new[[r]][cls$uid__])

  man <- if (length(manifest)) rbindlist(manifest) else
    data.table(rank = character(0), child = character(0), chosen_parent = character(0),
               n_parents = integer(0), n_records = integer(0), n_records_changed = integer(0),
               parents = character(0))
  list(classification = out[], manifest = man, n_rows_changed = n_rows_changed, n_orphans = n_orphans)
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
