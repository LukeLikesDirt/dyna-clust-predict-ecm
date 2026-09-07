#!/usr/bin/env Rscript
# remove_complexes.R — Detect and remove species-level "complexes" (species
# that cannot be discriminated by the barcode marker) before subsetting.
#
# This is NOT a parity port of dnabarcoder's -removecomplexes/removeComplexes.py.
# That algorithm (single-linkage cluster sequences at 100% identity; species
# co-occurring in a cluster form a complex; keep the species with most
# sequences, drop the rest) is unguarded against database contamination: one
# mislabeled/chimeric record with a spurious 100%-identity hit to an unrelated
# species is enough to transitively chain unrelated, well-separated species
# into one mega-complex. Measured on EUKARYOME full-ITS Aspergillus: unguarded
# single-linkage produced a 52-53 species complex spanning A. flavus, A.
# niger, A. fumigatus, A. nidulans -- species from different, well-separated
# sections with no real confusability, traced to one record whose 100%-
# identity partners span sections that share no real biology.
#
# This script adds one safeguard: before building the single-linkage complex
# graph, exclude any sequence whose 100%-identity matches span more than
# --max_hub_species distinct species. Validated on 6 fungal genera: this drops
# the largest complex from 52 -> 5 (Aspergillus), 39 -> 8 (Penicillium),
# 33 -> 17 (Trichoderma), 28 -> 8 (Alternaria), 11 -> 7 (Fusarium),
# 15 -> 6 (Cladosporium) -- sizes consistent with recognized real species
# complexes rather than contamination artifacts.
#
# Complexes are detected per-genus (a complex spans confusable species within
# one genus; cross-genus 100%-identity hits are themselves usually
# contamination, not biology, and are out of scope here).
#
# Deliberately NOT ported from dnabarcoder, per code review:
#   - MergeComplexes() reads a module-level `classes` global instead of taking
#     it as a parameter.
#   - Its `>=` tie-break on sequence count makes the kept representative
#     iteration-order-dependent. This script breaks ties on species name
#     (alphabetically first) for determinism.
#   - removeComplexes.py's own copy of ComputeFmeasure has a dead
#     `if (n==0) return f` guard before n is ever assigned.
#   - classname[-3]==" sp" in LoadClassification/GetTaxonName compares one
#     character to a 3-character string and is always False -- the intended
#     check (classname[-3:]==" sp") is not needed here since is_identified()
#     from R/utils.R already excludes "sp." placeholders correctly.
#
# Usage:
#   Rscript R/remove_complexes.R \
#     --fasta_in              data/full_ITS/eukaryome_ITS.fasta \
#     --classification_in     data/full_ITS/eukaryome_ITS.classification \
#     --fasta_out             data/full_ITS/eukaryome_ITS_nocomplex.fasta \
#     --classification_out    data/full_ITS/eukaryome_ITS_nocomplex.classification \
#     --manifest_out          data/full_ITS/complex_manifest.txt \
#     --threshold 1.0 --max_hub_species 3 --min_species_per_parent 2 \
#     --n_cpus 10 --tmp_dir ./tmp/remove_complexes
#
# Note: This script must be run from the project root directory, BEFORE
# 05_prepare_subsets.sh (subsetting should draw from complex-filtered data).

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(igraph)
  library(furrr)
  library(future)
})

source("R/utils.R")

# ── Arguments ─────────────────────────────────────────────────────────────────

option_list <- list(
  make_option("--fasta_in", type = "character", metavar = "FILE",
              help = "Input FASTA file [required]"),
  make_option("--classification_in", type = "character", metavar = "FILE",
              help = "Input tab-delimited classification file [required]"),
  make_option("--fasta_out", type = "character", metavar = "FILE",
              help = "Output FASTA with complex-member sequences removed [required]"),
  make_option("--classification_out", type = "character", metavar = "FILE",
              help = "Output classification, subset to surviving IDs [required]"),
  make_option("--manifest_out", type = "character",
              default = "complex_manifest.txt", metavar = "FILE",
              help = "Output path for the complex/hub manifest [default: %default]"),
  make_option("--rank", type = "character", default = "species", metavar = "STR",
              help = "Rank at which complexes are detected [default: %default]"),
  make_option("--parent_rank", type = "character", default = "genus", metavar = "STR",
              help = paste("Rank within which complexes are detected -- a complex is",
                           "a within-parent phenomenon (e.g. within-genus).",
                           "[default: %default]")),
  make_option("--threshold", type = "double", default = 1.0, metavar = "NUM",
              help = "Identity threshold defining a complex edge [default: %default]"),
  make_option("--max_hub_species", type = "integer", default = 3L, metavar = "INT",
              help = paste("Exclude a sequence from complex detection if its",
                           "threshold-identity matches span >= this many distinct",
                           "species -- guards against one mislabeled/contaminated",
                           "record chaining unrelated species together.",
                           "[default: %default]")),
  make_option("--min_species_per_parent", type = "integer", default = 2L, metavar = "INT",
              help = "Skip parent groups with fewer distinct species than this [default: %default]"),
  make_option("--require_complete", type = "character", default = "yes", metavar = "yes/no",
              help = paste("Restrict to rows with its_complete == TRUE (see",
                           "R/append_completeness.R) before complex detection, and drop",
                           "incomplete rows from the output entirely -- not just from",
                           "complex-detection evidence. Truncated/stub records inflate",
                           "similarity broadly (within-species too, not just cross-",
                           "species), so this is the gatekeeper for completeness as well",
                           "as complex removal. No-ops with a warning if the",
                           "classification file has no its_complete column.",
                           "[default: %default]")),
  make_option("--max_seq_no", type = "integer", default = 20000L, metavar = "INT",
              help = paste("Max sequences per parent group for complex detection;",
                           "excess is randomly downsampled (mirrors predict.R's own",
                           "cap) to bound per-group vsearch runtime/memory on the",
                           "largest genera. [default: %default]")),
  make_option("--iddef", type = "integer", default = 2L, metavar = "0-4",
              help = "vsearch --iddef pairwise identity definition [default: %default]"),
  make_option("--id_col", type = "character", default = "id", metavar = "STR",
              help = "ID column name in the classification file [default: %default]"),
  make_option("--n_cpus", type = "integer",
              default = max(1L, future::availableCores() - 1L), metavar = "INT",
              help = "Parallel workers [default: availableCores() - 1]"),
  make_option("--run_parallel", type = "character", default = "yes", metavar = "yes/no",
              help = "Process parent groups in parallel [default: %default]"),
  make_option("--tmp_dir", type = "character", default = "./tmp/remove_complexes",
              metavar = "DIR", help = "Directory for temporary vsearch output [default: %default]")
)

opt <- parse_args(
  OptionParser(
    option_list = option_list,
    usage = "%prog --fasta_in in.fa --classification_in in.tsv --fasta_out out.fa --classification_out out.tsv",
    description = paste(
      "Detect species-level complexes (indistinguishable by the barcode marker)",
      "via hub-guarded single-linkage clustering at --threshold, and remove all",
      "but the best-sampled species per complex."
    )
  )
)

required <- c("fasta_in", "classification_in", "fasta_out", "classification_out")
missing  <- required[sapply(required, function(x) is.null(opt[[x]]))]
if (length(missing) > 0) stop("Missing required options: ", paste0("--", missing, collapse = ", "))
if (!file.exists(opt$fasta_in))          stop("File not found: ", opt$fasta_in)
if (!file.exists(opt$classification_in)) stop("File not found: ", opt$classification_in)

rank                 <- opt$rank
parent_rank          <- opt$parent_rank
threshold            <- opt$threshold
max_hub_species      <- opt$max_hub_species
min_species_per_parent <- opt$min_species_per_parent
require_complete     <- identical(tolower(opt$require_complete), "yes")
max_seq_no           <- opt$max_seq_no
iddef                <- opt$iddef
id_col               <- opt$id_col
n_cpus               <- opt$n_cpus
run_parallel         <- identical(tolower(opt$run_parallel), "yes")
tmp_dir              <- opt$tmp_dir

for (d in unique(dirname(c(opt$fasta_out, opt$classification_out, opt$manifest_out)))) {
  if (!dir.exists(d) && d != ".") dir.create(d, recursive = TRUE)
}
if (!dir.exists(tmp_dir)) dir.create(tmp_dir, recursive = TRUE)

# ── Load inputs ───────────────────────────────────────────────────────────────

cat("[remove_complexes] Loading classification from:", opt$classification_in, "\n")
cls <- fread(opt$classification_in, sep = "\t", header = TRUE, colClasses = "character")
if (!id_col %in% names(cls))      stop("ID column '", id_col, "' not found.")
if (!rank %in% names(cls))        stop("Rank column '", rank, "' not found.")
if (!parent_rank %in% names(cls)) stop("Parent rank column '", parent_rank, "' not found.")

if (require_complete) {
  if ("its_complete" %in% names(cls)) {
    n_before <- nrow(cls)
    cls <- cls[tolower(its_complete) == "true"]
    cat(sprintf("[remove_complexes] --require_complete yes: restricted to its_complete rows: %d -> %d\n",
                n_before, nrow(cls)))
  } else {
    cat("[remove_complexes] WARNING: --require_complete yes but no its_complete column found",
        "(classification predates R/append_completeness.R) -- proceeding unfiltered.\n")
  }
}

cls_ok <- cls[is_identified(get(rank)) & is_identified(get(parent_rank))]
cat("[remove_complexes] Rows with identified", rank, "and", parent_rank, ":", nrow(cls_ok), "\n")

parent_groups <- split(cls_ok[[id_col]], cls_ok[[parent_rank]])
n_species_per_parent <- vapply(names(parent_groups), function(p) {
  uniqueN(cls_ok[[rank]][cls_ok[[id_col]] %in% parent_groups[[p]]])
}, integer(1))
parent_groups <- parent_groups[n_species_per_parent >= min_species_per_parent]
cat(sprintf("[remove_complexes] Parent groups with >= %d species: %d\n",
            min_species_per_parent, length(parent_groups)))

n_capped <- 0L
if (max_seq_no > 0) {
  parent_groups <- lapply(parent_groups, function(ids) {
    if (length(ids) > max_seq_no) {
      n_capped <<- n_capped + 1L
      sample(ids, max_seq_no)
    } else {
      ids
    }
  })
  if (n_capped > 0) {
    cat(sprintf("[remove_complexes] %d parent group(s) exceeded --max_seq_no %d, downsampled.\n",
                n_capped, max_seq_no))
  }
}

# ── Fixed-header FASTA index (for fast per-group subset extraction) ─────────

cat("[remove_complexes] Indexing FASTA:", opt$fasta_in, "\n")
fasta_lines <- readLines(opt$fasta_in)
h_idx       <- which(startsWith(fasta_lines, ">"))
h_ids       <- sub("^>([^ ]+).*", "\\1", fasta_lines[h_idx])
h_end       <- c(h_idx[-1] - 1L, length(fasta_lines))
id_to_pos   <- setNames(seq_along(h_idx), h_ids)

write_subset_fasta <- function(ids, out_path) {
  pos <- id_to_pos[ids]
  pos <- pos[!is.na(pos)]
  keep_lines <- unlist(lapply(pos, function(p) h_idx[p]:h_end[p]))
  writeLines(fasta_lines[keep_lines], out_path)
}

# ── Per-parent-group complex detection ────────────────────────────────────────
# Returns: list(removed_ids = character(), manifest_rows = data.table)

detect_complexes_one <- function(parent_name, ids, tmp_dir) {
  species_of <- setNames(cls_ok[[rank]][match(ids, cls_ok[[id_col]])], ids)
  seqcount   <- as.data.table(table(species_of))
  setnames(seqcount, c("species", "N"))

  safe_name <- gsub("[^A-Za-z0-9_]", "_", parent_name)
  fa_file   <- file.path(tmp_dir, paste0(safe_name, ".fasta"))
  sim_file  <- file.path(tmp_dir, paste0(safe_name, ".sim.txt"))
  write_subset_fasta(ids, fa_file)

  vsearch_cmd <- sprintf(
    "vsearch --allpairs_global '%s' --acceptall --userout '%s' --userfields query+target+id --iddef %d --threads 1",
    fa_file, sim_file, iddef
  )
  ok <- system(vsearch_cmd, ignore.stdout = TRUE, ignore.stderr = TRUE) == 0
  unlink(fa_file)
  if (!ok || !file.exists(sim_file) || file.info(sim_file)$size == 0) {
    unlink(sim_file)
    return(list(removed_ids = character(0), manifest_rows = data.table()))
  }

  sim <- fread(sim_file, header = FALSE, col.names = c("q", "t", "pid"),
              colClasses = c("character", "character", "numeric"))
  unlink(sim_file)
  sim[, score := pid / 100]
  sim <- sim[q %in% ids & t %in% ids]

  m <- sim[score >= threshold & q != t]
  m[, sq := species_of[q]][, st := species_of[t]]
  m <- m[sq != st]
  if (nrow(m) == 0) return(list(removed_ids = character(0), manifest_rows = data.table()))

  # ── Hub guard: flag records whose threshold-identity partners span too
  #    many distinct species, and exclude them from the complex graph ──────
  hub_score <- rbindlist(list(
    m[, .(id = q, partner_species = st)],
    m[, .(id = t, partner_species = sq)]
  ))[, .(n_distinct_partner_species = uniqueN(partner_species)), by = id]
  hub_ids <- hub_score[n_distinct_partner_species >= max_hub_species, id]

  m_guarded <- m[!(q %in% hub_ids) & !(t %in% hub_ids)]
  if (nrow(m_guarded) == 0) {
    hub_rows <- if (length(hub_ids) > 0) {
      hs <- hub_score[id %in% hub_ids]
      data.table(parent = parent_name, type = "hub_excluded", id = hs$id,
                species = species_of[hs$id], n_partner_species = hs$n_distinct_partner_species,
                complex_id = NA_character_, retained_representative = NA_character_)
    } else data.table()
    return(list(removed_ids = character(0), manifest_rows = hub_rows))
  }

  # ── Single-linkage complex graph over species names ──────────────────────
  edges <- unique(m_guarded[, .(a = pmin(sq, st), b = pmax(sq, st))])
  spp   <- sort(unique(species_of))
  g     <- graph_from_data_frame(edges, directed = FALSE, vertices = data.frame(name = spp))
  cm    <- components(g)

  dt <- data.table(species = spp, comp = cm$membership[spp])
  dt <- merge(dt, seqcount, by = "species", all.x = TRUE)
  dt[is.na(N), N := 0L]
  # Deterministic tie-break: most sequences wins; ties broken alphabetically
  # (dnabarcoder's own `>=` tie-break is iteration-order-dependent -- avoided here).
  setorder(dt, comp, -N, species)
  dt[, keep := seq_len(.N) == 1L, by = comp]

  complex_comps <- dt[, .N, by = comp][N > 1, comp]
  removed_species <- dt[comp %in% complex_comps & keep == FALSE, species]
  removed_ids <- names(species_of)[species_of %in% removed_species]

  manifest_rows <- rbindlist(list(
    if (nrow(dt[comp %in% complex_comps]) > 0) {
      cdt <- dt[comp %in% complex_comps]
      cdt[, `:=`(
        parent = parent_name, type = "complex_member",
        id = NA_character_, n_partner_species = NA_integer_,
        complex_id = paste0(safe_name, "_complex", sprintf("%02d", comp)),
        retained_representative = species[keep][1]
      ), by = comp][, .(parent, type, id, species, n_partner_species, complex_id, retained_representative)]
    } else data.table(),
    if (length(hub_ids) > 0) {
      hs <- hub_score[id %in% hub_ids]
      data.table(parent = parent_name, type = "hub_excluded", id = hs$id,
                species = species_of[hs$id], n_partner_species = hs$n_distinct_partner_species,
                complex_id = NA_character_, retained_representative = NA_character_)
    } else data.table()
  ), fill = TRUE)

  list(removed_ids = removed_ids, manifest_rows = manifest_rows)
}

# ── Run across all parent groups ──────────────────────────────────────────────

cat(sprintf("[remove_complexes] Detecting complexes within %d %s groups (threshold=%.4f, max_hub_species=%d)...\n",
            length(parent_groups), parent_rank, threshold, max_hub_species))

if (run_parallel) {
  if (.Platform$OS.type == "unix") plan(multicore, workers = n_cpus) else plan(multisession, workers = n_cpus)
} else {
  plan(sequential)
}

results <- future_map(
  names(parent_groups),
  function(p) detect_complexes_one(p, parent_groups[[p]], tmp_dir),
  .options = furrr_options(seed = TRUE),
  .progress = interactive()
)
plan(sequential)

all_removed_ids <- unlist(lapply(results, `[[`, "removed_ids"))
all_manifest    <- rbindlist(lapply(results, `[[`, "manifest_rows"), fill = TRUE)

cat(sprintf("[remove_complexes] %d sequences removed as complex members across %d parent groups.\n",
            length(all_removed_ids), length(parent_groups)))
cat(sprintf("[remove_complexes] %d records flagged and excluded as hubs.\n",
            sum(all_manifest$type == "hub_excluded", na.rm = TRUE)))

# ── Write outputs ─────────────────────────────────────────────────────────────

keep_ids <- setdiff(cls[[id_col]], all_removed_ids)
cls_out  <- cls[get(id_col) %in% keep_ids]
fwrite(cls_out, opt$classification_out, sep = "\t")
cat("[remove_complexes] Classification written to:", opt$classification_out,
    "(", nrow(cls_out), "rows )\n")

write_subset_fasta(keep_ids, opt$fasta_out)
cat("[remove_complexes] FASTA written to:", opt$fasta_out, "\n")

if (nrow(all_manifest) > 0) setorder(all_manifest, parent, type, species)
fwrite(all_manifest, opt$manifest_out, sep = "\t")
cat("[remove_complexes] Manifest written to:", opt$manifest_out, "\n")

cat("\n[remove_complexes] Done.\n")
