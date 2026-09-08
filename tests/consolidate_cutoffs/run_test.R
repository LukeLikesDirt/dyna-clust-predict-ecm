#!/usr/bin/env Rscript
# run_test.R — Correctness test for R/consolidate_cutoffs.R.
#
# Builds a small synthetic fixture (a raw nested cutoffs table + a matching
# classification file), runs the real consolidate_cutoffs.R script against it
# as a subprocess (exactly how the pipeline invokes it), then asserts the
# resolved output cell-by-cell against hand-traced expected values.
#
# Fixture design -- each scenario lives under its own, unrelated lineage root
# so the fallback chains never cross-contaminate:
#
#   KingdomA  control: full monotonic row, high confidence -> no change at all
#   KingdomB  gap at family; kingdom has no ancestor, so the gap can only be
#             filled from global
#   KingdomC  large monotonicity violation (order=0.90 > family=0.60);
#             confirms the clamp cascades (family AND genus both get raised)
#   KingdomD  small/noise violations (0.001-0.002 drops); confirms the clamp
#             is unconditional -- no tolerance exemption for small drops
#   KingdomE  species has a direct value but very low confidence; global
#             (lower cutoff, higher confidence) should win instead
#   OrderF/   ancestor-chain test: FamilyF's genus target is a gap with NO
#   FamilyF   direct kingdom/phylum/class rows above it either -- only its
#             immediate parent OrderF and global are candidates, and OrderF
#             should win on confidence. Also exercises a MIXED-source row
#             (OrderF's own family=global, genus=self) still getting clamped
#             correctly when genus < family after resolution.
#
# Usage: Rscript tests/consolidate_cutoffs/run_test.R
# Note:  Must be run from the project root directory.

suppressPackageStartupMessages(library(data.table))

fixture_dir <- "tests/consolidate_cutoffs/tmp"
dir.create(fixture_dir, showWarnings = FALSE, recursive = TRUE)

cutoffs_in  <- file.path(fixture_dir, "cutoffs.txt")
class_in    <- file.path(fixture_dir, "classification.txt")
output_path <- file.path(fixture_dir, "eukaryome_cutoffs.txt")

# ══════════════════════════════════════════════════════════════════════════════
# Build fixture: raw nested cutoffs table
# ══════════════════════════════════════════════════════════════════════════════

row <- function(rank, higher_rank, dataset, cutoff, confidence,
                seq_n = 100L, grp_n = 10L, max_prop = 0.2, multiseq_grp_n = NULL) {
  dt <- data.table(
    rank = rank, higher_rank = higher_rank, dataset = dataset,
    `cut-off` = cutoff, confidence = confidence,
    `sequence number` = seq_n, `group number` = grp_n, `max proportion` = max_prop
  )
  # Only rows in the KingdomG scenario set this explicitly, exercising the
  # singleton-aware override guard; every other scenario omits the column
  # entirely, exercising the multiseq_grp_n := grp_n backward-compat fallback
  # for cutoffs files produced before --min_multiseq_groups existed.
  if (!is.null(multiseq_grp_n)) dt[, `multiseq group number` := multiseq_grp_n]
  dt
}

cutoffs <- rbindlist(fill = TRUE, list(
  # ── Global anchor (monotonic, confidence 0.50 throughout) ──────────────────
  row("phylum",  "global", "All", 0.60, 0.50, 5000L, 50L,   0.10),
  row("class",   "global", "All", 0.68, 0.50, 5000L, 80L,   0.10),
  row("order",   "global", "All", 0.75, 0.50, 5000L, 150L,  0.10),
  row("family",  "global", "All", 0.80, 0.50, 5000L, 300L,  0.10),
  row("genus",   "global", "All", 0.85, 0.50, 5000L, 800L,  0.10),
  row("species", "global", "All", 0.95, 0.50, 5000L, 2000L, 0.10),

  # ── KingdomA: control, full monotonic row, high confidence ─────────────────
  row("phylum",  "kingdom", "KingdomA", 0.50, 0.90),
  row("class",   "kingdom", "KingdomA", 0.55, 0.90),
  row("order",   "kingdom", "KingdomA", 0.60, 0.90),
  row("family",  "kingdom", "KingdomA", 0.65, 0.90),
  row("genus",   "kingdom", "KingdomA", 0.70, 0.90),
  row("species", "kingdom", "KingdomA", 0.80, 0.90),

  # ── KingdomB: gap at family (no row at all) ─────────────────────────────────
  row("phylum",  "kingdom", "KingdomB", 0.55, 0.90),
  row("class",   "kingdom", "KingdomB", 0.60, 0.90),
  row("order",   "kingdom", "KingdomB", 0.65, 0.90),
  # family: deliberately absent
  row("genus",   "kingdom", "KingdomB", 0.83, 0.90),
  row("species", "kingdom", "KingdomB", 0.90, 0.90),

  # ── KingdomC: large monotonicity violation, self always wins on confidence ─
  row("phylum",  "kingdom", "KingdomC", 0.50, 0.90),
  row("class",   "kingdom", "KingdomC", 0.55, 0.90),
  row("order",   "kingdom", "KingdomC", 0.90, 0.90),
  row("family",  "kingdom", "KingdomC", 0.60, 0.90),
  row("genus",   "kingdom", "KingdomC", 0.85, 0.90),
  row("species", "kingdom", "KingdomC", 0.95, 0.90),

  # ── KingdomD: small/noise violations, self always wins on confidence ───────
  row("phylum",  "kingdom", "KingdomD", 0.860, 0.90),
  row("class",   "kingdom", "KingdomD", 0.859, 0.90),
  row("order",   "kingdom", "KingdomD", 0.858, 0.90),
  row("family",  "kingdom", "KingdomD", 0.870, 0.90),
  row("genus",   "kingdom", "KingdomD", 0.900, 0.90),
  row("species", "kingdom", "KingdomD", 0.950, 0.90),

  # ── KingdomE: species present but low confidence -> global should win ──────
  row("phylum",  "kingdom", "KingdomE", 0.50,  0.90),
  row("class",   "kingdom", "KingdomE", 0.55,  0.90),
  row("order",   "kingdom", "KingdomE", 0.60,  0.90),
  row("family",  "kingdom", "KingdomE", 0.65,  0.90),
  row("genus",   "kingdom", "KingdomE", 0.80,  0.90),
  row("species", "kingdom", "KingdomE", 0.999, 0.10),

  # ── OrderF / FamilyF: ancestor-chain test, no kingdom/phylum/class rows ────
  row("genus",   "order",  "OrderF",  0.77, 0.60, 400L, 40L, 0.15),
  row("species", "family", "FamilyF", 0.92, 0.85, 90L,  30L, 0.20),

  # ── GenusG1/FamilyG: singleton-aware override guard ─────────────────────────
  # self (GenusG1's own species prediction) has a modest but real confidence
  # from real, multi-sequence data (multiseq_grp_n=40 of 120 groups) -- above
  # the shared global anchor's 0.50 baseline so this scenario isolates the
  # guard's effect rather than being decided by global. Its ancestor (FamilyG's
  # species prediction) has HIGHER confidence but is almost entirely singletons
  # (multiseq_grp_n=8 of 60 groups) -- the classic singleton-inflation shape
  # (predict.R's --min_multiseq_groups guard keeps datasets this thin out of
  # production, but a dataset can still clear that bar with few multiseq
  # groups relative to a richer self row). Without the 2e guard the ancestor
  # wins purely on confidence; with it, the ancestor is excluded for bringing
  # less multi-sequence evidence than self.
  row("species", "genus",  "GenusG1", 0.85, 0.55, 300L, 120L, multiseq_grp_n = 40L),
  row("species", "family", "FamilyG", 0.99, 0.99, 60L,  60L,  multiseq_grp_n = 8L),

  # ── GenusH1/FamilyH: legitimate override still works ────────────────────────
  # Same shape, but the ancestor genuinely brings MORE multi-sequence evidence
  # than self (45 >= 40) as well as higher confidence -- confirms the guard
  # does not block a real, well-supported override.
  row("species", "genus",  "GenusH1", 0.85, 0.55, 300L, 120L, multiseq_grp_n = 40L),
  row("species", "family", "FamilyH", 0.99, 0.99, 500L, 500L, multiseq_grp_n = 45L)
))

fwrite(cutoffs, cutoffs_in, sep = "\t", quote = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
# Build fixture: classification file (for lineage lookup only)
# ══════════════════════════════════════════════════════════════════════════════

cls_row <- function(id, kingdom, phylum, class, order, family, genus, species) {
  data.table(id = id, kingdom = kingdom, phylum = phylum, class = class,
             order = order, family = family, genus = genus, species = species)
}

classification <- rbindlist(list(
  cls_row("s1", "KingdomA", "PhylumA", "ClassA", "OrderA", "FamilyA", "GenusA1", "SpeciesA1"),
  cls_row("s2", "KingdomB", "PhylumB", "ClassB", "OrderB", "FamilyB", "GenusB1", "SpeciesB1"),
  cls_row("s3", "KingdomC", "PhylumC", "ClassC", "OrderC", "FamilyC", "GenusC1", "SpeciesC1"),
  cls_row("s4", "KingdomD", "PhylumD", "ClassD", "OrderD", "FamilyD", "GenusD1", "SpeciesD1"),
  cls_row("s5", "KingdomE", "PhylumE", "ClassE", "OrderE", "FamilyE", "GenusE1", "SpeciesE1"),
  # F-lineage: 3 consistent rows so the modal parent is unambiguous
  cls_row("f1", "KingdomF", "PhylumF", "ClassF", "OrderF", "FamilyF", "GenusF1", "SpeciesF1"),
  cls_row("f2", "KingdomF", "PhylumF", "ClassF", "OrderF", "FamilyF", "GenusF2", "SpeciesF2"),
  cls_row("f3", "KingdomF", "PhylumF", "ClassF", "OrderF", "FamilyF", "GenusF1", "SpeciesF3"),
  cls_row("g1", "KingdomG", "PhylumG", "ClassG", "OrderG", "FamilyG", "GenusG1", "SpeciesG1"),
  cls_row("h1", "KingdomH", "PhylumH", "ClassH", "OrderH", "FamilyH", "GenusH1", "SpeciesH1")
))

fwrite(classification, class_in, sep = "\t", quote = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
# Run the real script as a subprocess
# ══════════════════════════════════════════════════════════════════════════════

cat("Running R/consolidate_cutoffs.R on the synthetic fixture...\n\n")

rc <- system2(
  "Rscript",
  args = c("R/consolidate_cutoffs.R",
           "--cutoffs_in", cutoffs_in,
           "--classification_in", class_in,
           "--output", output_path),
  stdout = TRUE, stderr = TRUE
)
cat(paste(rc, collapse = "\n"), "\n\n")

exit_status <- attr(rc, "status")
if (!is.null(exit_status) && exit_status != 0) {
  cat("*** FAIL: consolidate_cutoffs.R exited with status", exit_status, "***\n")
  quit(status = 1)
}

if (!file.exists(output_path)) {
  cat("*** FAIL: output file was not created:", output_path, "***\n")
  quit(status = 1)
}

result <- fread(output_path, sep = "\t", header = TRUE)
setnames(result,
         old = c("cut-off", "sequence number", "group number", "max proportion"),
         new = c("cutoff", "seq_n", "grp_n", "max_prop"))

# ══════════════════════════════════════════════════════════════════════════════
# Assertions
# ══════════════════════════════════════════════════════════════════════════════

n_pass <- 0L
n_fail <- 0L

get_cell <- function(hr, ds_name, target_rank) {
  # Parameter named ds_name (not "dataset") to avoid colliding with the
  # `dataset` column inside data.table's [ NSE.
  hit <- result[higher_rank == hr & dataset == ds_name & rank == target_rank]
  if (nrow(hit) == 0) return(NULL)
  hit[1]
}

check <- function(label, hr, ds_name, target_rank, field, expected, tol = 1e-9) {
  cell <- get_cell(hr, ds_name, target_rank)
  if (is.null(cell)) {
    cat(sprintf("FAIL  %-55s  row not found (%s / %s / %s)\n", label, hr, ds_name, target_rank))
    n_fail <<- n_fail + 1L
    return(invisible())
  }
  actual <- cell[[field]]
  ok <- if (is.na(expected)) {
    is.na(actual)
  } else if (is.numeric(expected)) {
    !is.na(actual) && abs(actual - expected) < tol
  } else {
    identical(as.character(actual), as.character(expected))
  }
  if (ok) {
    cat(sprintf("PASS  %-55s  %s = %s\n", label, field, actual))
    n_pass <<- n_pass + 1L
  } else {
    cat(sprintf("FAIL  %-55s  %s = %s (expected %s)\n", label, field, actual, expected))
    n_fail <<- n_fail + 1L
  }
}

cat("══════════════════════════════════════════════════════════════\n")
cat("KingdomA -- control: no gaps, no violations, self always wins\n")
cat("══════════════════════════════════════════════════════════════\n")
for (r in list(c("phylum", 0.50), c("class", 0.55), c("order", 0.60),
               c("family", 0.65), c("genus", 0.70), c("species", 0.80))) {
  check(sprintf("KingdomA.%s cutoff unchanged", r[1]), "kingdom", "KingdomA", r[1], "cutoff", as.numeric(r[2]))
  check(sprintf("KingdomA.%s source is self", r[1]), "kingdom", "KingdomA", r[1], "source", "self")
  check(sprintf("KingdomA.%s not clamped", r[1]), "kingdom", "KingdomA", r[1], "clamped", "FALSE")
}

cat("\n══════════════════════════════════════════════════════════════\n")
cat("KingdomB -- gap at family, filled from global (kingdom has no ancestor)\n")
cat("══════════════════════════════════════════════════════════════\n")
check("KingdomB.family cutoff from global", "kingdom", "KingdomB", "family", "cutoff", 0.80)
check("KingdomB.family source is global", "kingdom", "KingdomB", "family", "source", "global")
check("KingdomB.family confidence from global", "kingdom", "KingdomB", "family", "confidence", 0.50)
check("KingdomB.family original_cutoff is NA (true gap)", "kingdom", "KingdomB", "family", "original_cutoff", NA_real_)
check("KingdomB.family not clamped", "kingdom", "KingdomB", "family", "clamped", "FALSE")
check("KingdomB.genus unaffected", "kingdom", "KingdomB", "genus", "cutoff", 0.83)
check("KingdomB.species unaffected", "kingdom", "KingdomB", "species", "cutoff", 0.90)

cat("\n══════════════════════════════════════════════════════════════\n")
cat("KingdomC -- large violation (order=0.90 > family=0.60): cascading clamp\n")
cat("══════════════════════════════════════════════════════════════\n")
check("KingdomC.phylum unclamped", "kingdom", "KingdomC", "phylum", "cutoff", 0.50)
check("KingdomC.class unclamped", "kingdom", "KingdomC", "class", "cutoff", 0.55)
check("KingdomC.order unclamped", "kingdom", "KingdomC", "order", "cutoff", 0.90)
check("KingdomC.family clamped up to order's 0.90", "kingdom", "KingdomC", "family", "cutoff", 0.90)
check("KingdomC.family clamped flag is TRUE", "kingdom", "KingdomC", "family", "clamped", "TRUE")
check("KingdomC.family original_cutoff preserved", "kingdom", "KingdomC", "family", "original_cutoff", 0.60)
check("KingdomC.genus ALSO clamped up to 0.90 (cascade)", "kingdom", "KingdomC", "genus", "cutoff", 0.90)
check("KingdomC.genus clamped flag is TRUE", "kingdom", "KingdomC", "genus", "clamped", "TRUE")
check("KingdomC.genus original_cutoff preserved", "kingdom", "KingdomC", "genus", "original_cutoff", 0.85)
check("KingdomC.species unclamped", "kingdom", "KingdomC", "species", "cutoff", 0.95)

cat("\n══════════════════════════════════════════════════════════════\n")
cat("KingdomD -- small/noise violations still get clamped (no tolerance)\n")
cat("══════════════════════════════════════════════════════════════\n")
check("KingdomD.phylum unclamped", "kingdom", "KingdomD", "phylum", "cutoff", 0.860)
check("KingdomD.class clamped up to 0.860", "kingdom", "KingdomD", "class", "cutoff", 0.860)
check("KingdomD.class clamped flag is TRUE", "kingdom", "KingdomD", "class", "clamped", "TRUE")
check("KingdomD.order clamped up to 0.860", "kingdom", "KingdomD", "order", "cutoff", 0.860)
check("KingdomD.order clamped flag is TRUE", "kingdom", "KingdomD", "order", "clamped", "TRUE")
check("KingdomD.family unclamped", "kingdom", "KingdomD", "family", "cutoff", 0.870)
check("KingdomD.genus unclamped", "kingdom", "KingdomD", "genus", "cutoff", 0.900)
check("KingdomD.species unclamped", "kingdom", "KingdomD", "species", "cutoff", 0.950)

cat("\n══════════════════════════════════════════════════════════════\n")
cat("KingdomE -- low-confidence direct value loses to higher-confidence global\n")
cat("══════════════════════════════════════════════════════════════\n")
check("KingdomE.genus unaffected (self wins, high confidence)", "kingdom", "KingdomE", "genus", "cutoff", 0.80)
check("KingdomE.species overridden by global", "kingdom", "KingdomE", "species", "cutoff", 0.95)
check("KingdomE.species source is global", "kingdom", "KingdomE", "species", "source", "global")
check("KingdomE.species original_cutoff preserved", "kingdom", "KingdomE", "species", "original_cutoff", 0.999)
check("KingdomE.species original_confidence preserved", "kingdom", "KingdomE", "species", "original_confidence", 0.10)
check("KingdomE.species not clamped (0.95 already >= 0.80)", "kingdom", "KingdomE", "species", "clamped", "FALSE")

cat("\n══════════════════════════════════════════════════════════════\n")
cat("FamilyF -- ancestor-chain fallback: OrderF beats global on confidence\n")
cat("══════════════════════════════════════════════════════════════\n")
check("FamilyF.genus resolved from parent OrderF", "family", "FamilyF", "genus", "cutoff", 0.77)
check("FamilyF.genus source is order:OrderF", "family", "FamilyF", "genus", "source", "order:OrderF")
check("FamilyF.genus original_cutoff is NA (true gap)", "family", "FamilyF", "genus", "original_cutoff", NA_real_)
check("FamilyF.species self wins over global", "family", "FamilyF", "species", "cutoff", 0.92)
check("FamilyF.species source is self", "family", "FamilyF", "species", "source", "self")

cat("\n══════════════════════════════════════════════════════════════\n")
cat("OrderF -- mixed-source row (family=global, genus=self) still clamps\n")
cat("══════════════════════════════════════════════════════════════\n")
check("OrderF.family filled from global (no ancestor above OrderF)", "order", "OrderF", "family", "cutoff", 0.80)
check("OrderF.family source is global", "order", "OrderF", "family", "source", "global")
check("OrderF.genus clamped up from 0.77 to family's 0.80", "order", "OrderF", "genus", "cutoff", 0.80)
check("OrderF.genus clamped flag is TRUE", "order", "OrderF", "genus", "clamped", "TRUE")
check("OrderF.genus source remains self (clamp doesn't change source)", "order", "OrderF", "genus", "source", "self")
check("OrderF.genus original_cutoff preserved", "order", "OrderF", "genus", "original_cutoff", 0.77)
check("OrderF.species filled from global", "order", "OrderF", "species", "cutoff", 0.95)

cat("\n══════════════════════════════════════════════════════════════\n")
cat("GenusG1 -- singleton-aware guard: weak-evidence ancestor excluded\n")
cat("══════════════════════════════════════════════════════════════\n")
check("GenusG1.species self wins despite lower confidence", "genus", "GenusG1", "species", "cutoff", 0.85)
check("GenusG1.species source is self (ancestor excluded by guard)", "genus", "GenusG1", "species", "source", "self")
check("GenusG1.species confidence is self's own (0.55), not ancestor's 0.99", "genus", "GenusG1", "species", "confidence", 0.55)

cat("\n══════════════════════════════════════════════════════════════\n")
cat("GenusH1 -- guard does not block a genuinely well-supported override\n")
cat("══════════════════════════════════════════════════════════════\n")
check("GenusH1.species overridden by ancestor (>= self's multiseq evidence)", "genus", "GenusH1", "species", "cutoff", 0.99)
check("GenusH1.species source is family:FamilyH", "genus", "GenusH1", "species", "source", "family:FamilyH")
check("GenusH1.species confidence is ancestor's 0.99", "genus", "GenusH1", "species", "confidence", 0.99)

cat("\n══════════════════════════════════════════════════════════════\n")
cat(sprintf("RESULT: %d passed, %d failed\n", n_pass, n_fail))
cat("══════════════════════════════════════════════════════════════\n")

if (n_fail > 0) quit(status = 1)
