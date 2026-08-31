#!/usr/bin/env Rscript
# consolidate_cutoffs.R — Fill gaps and repair monotonicity in a region's
# nested similarity-cutoff table.
#
# Two problems in the raw <prefix>.cutoffs.json.txt produced by predict.R
# (via 06a_predict_cutoffs.sh / 06b_predict_cutoffs_region.sh):
#
#   1. Gaps: a (higher_rank, dataset, rank) cell is missing whenever that
#      dataset's subset failed subset.R's min_subgroups/min_sequences filters,
#      so predict.R never attempted a computation for it.
#   2. Monotonicity violations: sequence identity must increase from
#      phylum -> species (each rank nests inside the one above it), but two
#      independently-computed ranks can disagree.
#
# For every (higher_rank, dataset) that has at least one direct computation,
# and every target rank valid under that higher_rank, this script:
#
#   1. Builds a candidate list: the dataset's own direct value (if any), each
#      ancestor's direct value at the same target rank (walking the real
#      taxonomic lineage up to kingdom, derived from the classification
#      file), and the eukaryome-wide global value for that target rank.
#      Only directly-observed rows are used at each level -- never an
#      already-resolved value from this same pass -- so there is no
#      fallback-of-a-fallback circularity.
#   2. Picks the candidate with the highest confidence (F-measure). This one
#      rule both fills gaps (only ancestors/global exist) and overrides a
#      weak direct value (an ancestor or global beats it on confidence).
#      Ties break toward the more specific (deeper) candidate.
#   3. Clamps each dataset's own resolved row to be non-decreasing from its
#      coarsest to its finest target rank (cumulative max), unconditionally
#      -- this can raise a value that was already the "winner" in step 2.
#
# Usage:
#   Rscript R/consolidate_cutoffs.R \
#     --cutoffs_in        data/full_ITS/eukaryome.cutoffs.json.txt \
#     --classification_in data/full_ITS/eukaryome_ITS.classification \
#     --output            data/full_ITS/eukaryome_cutoffs.txt
#
# Note: This script must be run from the project root directory.

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
})

source("R/utils.R")

# ── Arguments ─────────────────────────────────────────────────────────────────

option_list <- list(
  make_option("--cutoffs_in",
              type = "character", metavar = "FILE",
              help = "Raw <prefix>.cutoffs.json.txt for one region [required]"),
  make_option("--classification_in",
              type = "character", metavar = "FILE",
              help = "Region classification file, for lineage lookup only [required]"),
  make_option("--output",
              type = "character", metavar = "FILE",
              help = "Output path for the consolidated table [required]")
)

opt <- parse_args(
  OptionParser(
    option_list = option_list,
    usage       = "%prog --cutoffs_in FILE --classification_in FILE --output FILE",
    description = paste(
      "Fill gaps and repair monotonicity in a region's nested cutoff table",
      "via a confidence-ranked fallback chain (self -> ancestors -> global)."
    )
  )
)

if (is.null(opt$cutoffs_in))        stop("--cutoffs_in is required.")
if (is.null(opt$classification_in)) stop("--classification_in is required.")
if (is.null(opt$output))            stop("--output is required.")
if (!file.exists(opt$cutoffs_in))        stop("File not found: ", opt$cutoffs_in)
if (!file.exists(opt$classification_in)) stop("File not found: ", opt$classification_in)

GLOBAL_HIGHER_RANK <- "global"
GLOBAL_DATASET     <- "All"   # matches predict.R's generate_datasets(): list(All = ids)

# Ranks that ever act as a parent (kingdom..genus); species is never a parent.
PARENT_RANKS <- rank_hierarchy[rank_hierarchy != "species"]

# ── Load inputs ───────────────────────────────────────────────────────────────

cat("[consolidate] Loading cutoffs from:", opt$cutoffs_in, "\n")
cutoffs <- fread(opt$cutoffs_in, sep = "\t", header = TRUE)
setnames(cutoffs,
         old = c("cut-off", "sequence number", "group number", "max proportion"),
         new = c("cutoff", "seq_n", "grp_n", "max_prop"))

cat("[consolidate] Loading classification from:", opt$classification_in, "\n")
# Default quoting (not quote = "") -- see R/predict.R load_classification()
# for why this matters for data written by dereplicate_lca.R's fwrite().
cls <- fread(opt$classification_in, sep = "\t", header = TRUE,
             fill = TRUE, na.strings = "", data.table = TRUE)

if (nrow(cutoffs[higher_rank == GLOBAL_HIGHER_RANK]) == 0) {
  warning("No '", GLOBAL_HIGHER_RANK, "' rows found in ", opt$cutoffs_in,
          " -- the fallback chain has no top-level anchor. Run the global ",
          "prediction step first (06a/06b, after the *_pred_id_global.txt ",
          "filename fix).", call. = FALSE)
}

# ── Build the taxonomic lineage table ────────────────────────────────────────
# For each rank phylum..genus, map each taxon name to its parent at the rank
# above, taken as the mode across the classification file (handles the rare
# case of an inconsistent LCA assignment for one name).

build_lineage_table <- function(cls) {
  lineage <- list()
  child_ranks <- rank_hierarchy[2:6]   # phylum, class, order, family, genus
  for (child_rank in child_ranks) {
    parent_rank <- rank_hierarchy[match(child_rank, rank_hierarchy) - 1]
    sub <- cls[is_identified(get(child_rank)) & is_identified(get(parent_rank)),
               .(child = get(child_rank), parent = get(parent_rank))]
    if (nrow(sub) == 0) { lineage[[child_rank]] <- character(0); next }

    counts <- sub[, .N, by = .(child, parent)]
    setorder(counts, child, -N)
    modal <- counts[, .SD[1], by = child]

    total_by_child <- counts[, .(total = sum(N)), by = child]
    check <- merge(modal, total_by_child, by = "child")
    inconsistent <- check[N / total < 0.9]
    if (nrow(inconsistent) > 0) {
      warning(sprintf(
        paste("[consolidate] %d %s name(s) have an inconsistent parent %s across",
              "the classification file (modal parent used for lineage lookup,",
              "e.g. %s)."),
        nrow(inconsistent), child_rank, parent_rank, inconsistent$child[1]
      ), call. = FALSE)
    }

    lineage[[child_rank]] <- setNames(modal$parent, modal$child)
  }
  lineage
}

cat("[consolidate] Building taxonomic lineage table...\n")
lineage <- build_lineage_table(cls)

# ── Ancestor chain for one (rank, taxon) pair, nearest first, up to kingdom ──

get_ancestor_chain <- function(rank, taxon, lineage) {
  chain <- list()
  current_rank <- rank
  current_name <- taxon
  while (current_rank != "kingdom") {
    parent_rank <- rank_hierarchy[match(current_rank, rank_hierarchy) - 1]
    parent_name <- lineage[[current_rank]][[current_name]]
    if (is.null(parent_name) || is.na(parent_name)) break
    chain[[length(chain) + 1]] <- list(rank = parent_rank, name = parent_name)
    current_rank <- parent_rank
    current_name <- parent_name
  }
  chain
}

# ── Fast lookup: (higher_rank, dataset, rank) -> row ─────────────────────────

KEY_SEP <- "\x1f"  # unit separator; unambiguous regardless of taxon name content
cutoffs[, lookup_key := paste(higher_rank, dataset, rank, sep = KEY_SEP)]
cutoffs_lookup <- split(cutoffs, seq_len(nrow(cutoffs)))
names(cutoffs_lookup) <- cutoffs$lookup_key

lookup_row <- function(hr, ds, rk) {
  key <- paste(hr, ds, rk, sep = KEY_SEP)
  cutoffs_lookup[[key]]
}

# ── Resolve one (higher_rank, dataset, target_rank) cell ─────────────────────
# Returns a list describing the winning candidate, or NULL if nothing at all
# (not even global) was found -- reported, not silently dropped.

resolve_cell <- function(hr, dataset, target_rank, lineage) {
  candidates <- list()

  self_row <- lookup_row(hr, dataset, target_rank)
  if (!is.null(self_row)) {
    candidates[[length(candidates) + 1]] <- list(
      source = "self", cutoff = self_row$cutoff, confidence = self_row$confidence,
      seq_n = self_row$seq_n, grp_n = self_row$grp_n, max_prop = self_row$max_prop
    )
  }

  for (anc in get_ancestor_chain(hr, dataset, lineage)) {
    anc_row <- lookup_row(anc$rank, anc$name, target_rank)
    if (!is.null(anc_row)) {
      candidates[[length(candidates) + 1]] <- list(
        source = sprintf("%s:%s", anc$rank, anc$name),
        cutoff = anc_row$cutoff, confidence = anc_row$confidence,
        seq_n = anc_row$seq_n, grp_n = anc_row$grp_n, max_prop = anc_row$max_prop
      )
    }
  }

  global_row <- lookup_row(GLOBAL_HIGHER_RANK, GLOBAL_DATASET, target_rank)
  if (!is.null(global_row)) {
    candidates[[length(candidates) + 1]] <- list(
      source = "global", cutoff = global_row$cutoff, confidence = global_row$confidence,
      seq_n = global_row$seq_n, grp_n = global_row$grp_n, max_prop = global_row$max_prop
    )
  }

  candidates <- Filter(function(c) !is.na(c$confidence), candidates)
  if (length(candidates) == 0) return(NULL)

  # which.max returns the FIRST maximum, and candidates are ordered
  # self -> nearest ancestor -> ... -> global, so ties break toward
  # specificity for free.
  best <- candidates[[which.max(vapply(candidates, `[[`, numeric(1), "confidence"))]]

  list(
    source            = best$source,
    cutoff            = best$cutoff,
    confidence        = best$confidence,
    seq_n             = best$seq_n,
    grp_n             = best$grp_n,
    max_prop          = best$max_prop,
    original_cutoff   = if (!is.null(self_row)) self_row$cutoff     else NA_real_,
    original_confidence = if (!is.null(self_row)) self_row$confidence else NA_real_
  )
}

# ── Enumerate cells and resolve ──────────────────────────────────────────────

cat("[consolidate] Resolving cells by confidence-ranked fallback chain...\n")

results <- list()
n_unresolved <- 0L

for (hr in PARENT_RANKS) {
  target_ranks <- rank_hierarchy[(match(hr, rank_hierarchy) + 1):length(rank_hierarchy)]
  datasets_at_hr <- unique(cutoffs[higher_rank == hr]$dataset)
  if (length(datasets_at_hr) == 0) next

  for (ds in datasets_at_hr) {
    for (target_rank in target_ranks) {
      resolved <- resolve_cell(hr, ds, target_rank, lineage)
      if (is.null(resolved)) {
        n_unresolved <- n_unresolved + 1L
        results[[length(results) + 1]] <- data.table(
          rank = target_rank, higher_rank = hr, dataset = ds,
          cutoff = NA_real_, confidence = NA_real_,
          seq_n = NA_integer_, grp_n = NA_integer_, max_prop = NA_real_,
          source = "unresolved", clamped = FALSE,
          original_cutoff = NA_real_, original_confidence = NA_real_
        )
        next
      }
      results[[length(results) + 1]] <- data.table(
        rank = target_rank, higher_rank = hr, dataset = ds,
        cutoff = resolved$cutoff, confidence = resolved$confidence,
        seq_n = resolved$seq_n, grp_n = resolved$grp_n, max_prop = resolved$max_prop,
        source = resolved$source, clamped = FALSE,
        original_cutoff = resolved$original_cutoff,
        original_confidence = resolved$original_confidence
      )
    }
  }
}

consolidated <- rbindlist(results)

if (n_unresolved > 0) {
  cat(sprintf("[consolidate] WARNING: %d cell(s) had no candidate at any level (not even global).\n",
              n_unresolved))
}

# ── Monotonicity clamp: cumulative max, per (higher_rank, dataset) row ───────
# Applied unconditionally regardless of violation size (no tolerance check).

cat("[consolidate] Applying cumulative-max monotonicity clamp...\n")

consolidated[, rank := factor(rank, levels = rank_hierarchy)]
setorder(consolidated, higher_rank, dataset, rank)

clamp_group <- function(cutoff) {
  out <- cutoff
  running_max <- -Inf
  for (i in seq_along(out)) {
    if (is.na(out[i])) next
    if (out[i] < running_max) out[i] <- running_max
    running_max <- max(running_max, out[i])
  }
  out
}

consolidated[, clamped_cutoff := clamp_group(cutoff), by = .(higher_rank, dataset)]
consolidated[, clamped := !is.na(cutoff) & !is.na(clamped_cutoff) & clamped_cutoff != cutoff]
consolidated[clamped == TRUE, cutoff := clamped_cutoff]
consolidated[, clamped_cutoff := NULL]

n_clamped <- sum(consolidated$clamped)
cat(sprintf("[consolidate] %d cell(s) raised by the monotonicity clamp.\n", n_clamped))

# ── Save ──────────────────────────────────────────────────────────────────────

setnames(consolidated,
         old = c("cutoff", "seq_n", "grp_n", "max_prop"),
         new = c("cut-off", "sequence number", "group number", "max proportion"))
setcolorder(consolidated, c(
  "rank", "higher_rank", "dataset", "cut-off", "confidence",
  "sequence number", "group number", "max proportion",
  "source", "clamped", "original_cutoff", "original_confidence"
))

fwrite(consolidated, opt$output, sep = "\t", quote = FALSE, na = "")

cat("\n[consolidate] Summary by source:\n")
print(consolidated[, .N, by = source][order(-N)])

cat("\n[consolidate] Consolidated table written to:", opt$output, "\n")
cat(sprintf("[consolidate]   %d rows | %d clamped | %d unresolved\n",
            nrow(consolidated), n_clamped, n_unresolved))
