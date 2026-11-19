#!/usr/bin/env Rscript
# test_iddef.R — Regression test for the identity definition (--iddef, default
# 2; --iddef 1 on request) and --redo yes in R/predict.R and R/compute_sim.R.
# The minimum overlap that goes with --iddef 2 (--min_overlap, in bp) is tested
# in tests/test_min_overlap.R.
# Replaces test_min_cov.R: the --min_cov coverage gate was removed on
# 2026-11-10, when the pipeline moved from --iddef 2 + gate to --iddef 1, and
# the pipeline went back to --iddef 2 on 2026-11-18 with a bp overlap in place
# of a fractional gate. Both definitions stay selectable, so both are tested.
#
# The fixture is the case the gate existed for. Every sequence of species A
# ends with the 30 bases that every sequence of species B starts with. --iddef 2
# excludes terminal gaps, so it scores each A-B pair 100% from those 30 bases;
# --iddef 1 counts the terminal gaps, so the same pair scores about 15%.
# Within a species, sequences differ by one substitution from a common parent
# (at least 99% identical to each other).
#
# So at threshold 0.99 the correct clustering is two clusters and F = 1, while a
# spurious A-B link merges them and gives F = 2 x 4 / (4 + 8) = 0.6667.
#
# Runs predict.R as a subprocess, because sourcing it would trigger its CLI
# parsing.
#
# Usage: Rscript tests/test_iddef.R
# Note:  Must be run from the project root directory, with vsearch on PATH.

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

root <- tempfile("test_iddef_")
dir.create(root)

n_pass <- 0L
n_fail <- 0L
check <- function(desc, ok) {
  if (isTRUE(ok)) {
    n_pass <<- n_pass + 1L
    cat("PASS ", desc, "\n", sep = "")
  } else {
    n_fail <<- n_fail + 1L
    cat("FAIL ", desc, "\n", sep = "")
  }
}

# ── Fixture ───────────────────────────────────────────────────────────────────

set.seed(42)
bases      <- c("A", "C", "G", "T")
random_seq <- function(n) paste(sample(bases, n, replace = TRUE), collapse = "")

# Parent plus three variants, each with one substitution at its own position.
species_group <- function(parent, name) {
  chars    <- strsplit(parent, "")[[1]]
  variants <- vapply(c(60L, 100L, 140L), function(pos) {
    v <- chars
    v[pos] <- setdiff(bases, v[pos])[1]
    paste(v, collapse = "")
  }, character(1))
  setNames(c(parent, variants), paste0(name, "_", 1:4))
}

shared <- random_seq(30)
A <- species_group(paste0(random_seq(170), shared), "a")
B <- species_group(paste0(shared, random_seq(170)), "b")

write_fixture <- function(seqs, species, stem) {
  fasta <- paste0(stem, ".fasta")
  cls   <- paste0(stem, ".classification")
  writeLines(paste0(">", names(seqs), "\n", seqs), fasta)
  fwrite(data.table(id = names(seqs), genus = "Testgenus", species = species), cls, sep = "\t")
  list(fasta = fasta, cls = cls)
}

fx_ab <- write_fixture(c(A, B), rep(c("Testgenus_alpha", "Testgenus_beta"), each = 4),
                       file.path(root, "ab"))

# ── Helpers ───────────────────────────────────────────────────────────────────

run_predict <- function(fx, out, ..., higher_rank = NULL) {
  dir.create(file.path(out, "tmp"), showWarnings = FALSE)
  args <- c("R/predict.R", "--input", fx$fasta, "--classification", fx$cls,
            "--rank", "species", "--start_threshold", "0.9", "--end_threshold", "1",
            "--step", "0.01", "--prefix", "t", "--out", out,
            "--tmp_dir", file.path(out, "tmp"),
            "--run_parallel", "no", "--n_cpus", "1",
            if (!is.null(higher_rank)) c("--higher_rank", higher_rank), ...)
  log    <- suppressWarnings(system2("Rscript", args, stdout = TRUE, stderr = TRUE))
  status <- attr(log, "status")
  list(status = if (is.null(status)) 0L else status, log = log)
}

# F-measure recorded in <out>/t.predicted for one dataset and threshold.
f_at <- function(out, threshold = "0.9900", dataset = "All") {
  pred <- fromJSON(file.path(out, "t.predicted"), simplifyVector = FALSE)
  as.numeric(pred$species[[dataset]]$fmeasures[[threshold]])
}

new_dir <- function(name) {
  d <- file.path(root, name)
  dir.create(d)
  d
}

all_pairs <- function(fasta, iddef) {
  hits <- file.path(root, paste0("hits_iddef", iddef, ".txt"))
  system2("vsearch", c("--allpairs_global", fasta, "--acceptall", "--userout", hits,
                       "--userfields", "query+target+id+qcov+tcov", "--iddef", iddef,
                       "--threads", "1"), stdout = FALSE, stderr = FALSE)
  fread(hits, header = FALSE, col.names = c("q", "t", "id", "qcov", "tcov"))
}

# ── The fixture must actually trigger the artefact, or the rest proves nothing ─

cross2 <- all_pairs(fx_ab$fasta, 2)[substr(q, 1, 1) != substr(t, 1, 1)]
check("fixture: under --iddef 2 every A-B pair is a 100% match on under 20% coverage",
      nrow(cross2) == 16L && all(cross2$id == 100) && all(cross2$qcov < 20) && all(cross2$tcov < 20))
cross1 <- all_pairs(fx_ab$fasta, 1)[substr(q, 1, 1) != substr(t, 1, 1)]
check("fixture: under --iddef 1 the same 16 pairs score under 20%",
      nrow(cross1) == 16L && all(cross1$id < 20))

# ── --iddef 1: global prediction ──────────────────────────────────────────────

d_def <- new_dir("iddef1")
r_def <- run_predict(fx_ab, d_def, "--iddef", "1")
check("--iddef 1 run succeeds", r_def$status == 0L)
check("--iddef 1 (no coverage filter): no spurious link, F = 1", f_at(d_def) == 1)

d_two_def <- new_dir("iddef2")
r_two_def <- run_predict(fx_ab, d_two_def, "--iddef", "2")
check("--iddef 2 run succeeds", r_two_def$status == 0L)
check("--iddef 2: the 100% artefact merges the species, F = 0.6667", f_at(d_two_def) == 0.6667)

d_dflt <- new_dir("default")
invisible(run_predict(fx_ab, d_dflt))
check("the default is --iddef 2 (same F = 0.6667 as an explicit --iddef 2)", f_at(d_dflt) == 0.6667)

# ── Local prediction, which goes through process_dataset() and workers ────────

d_local <- new_dir("local")
r_local <- run_predict(fx_ab, d_local, "--iddef", "1", higher_rank = "genus")
check("local run succeeds", r_local$status == 0L)
check("local: per-dataset vsearch uses --iddef 1, F = 1",
      f_at(d_local, dataset = "Testgenus") == 1)

d_par <- new_dir("parallel")
r_par <- run_predict(fx_ab, d_par, "--iddef", "1", higher_rank = "genus",
                     "--run_parallel", "yes", "--n_cpus", "2")
check("parallel local run succeeds", r_par$status == 0L)
check("parallel: the workers use --iddef 1, F = 1", f_at(d_par, dataset = "Testgenus") == 1)

d_par2 <- new_dir("parallel_iddef2")
r_par2 <- run_predict(fx_ab, d_par2, higher_rank = "genus", "--run_parallel", "yes", "--n_cpus", "2",
                      "--iddef", "2")
check("parallel local run with --iddef 2 succeeds", r_par2$status == 0L)
check("parallel: --iddef reaches the workers (--iddef 2 gives F = 0.6667)",
      f_at(d_par2, dataset = "Testgenus") == 0.6667)

# ── The coverage gate is gone ─────────────────────────────────────────────────

d_gone <- new_dir("gone")
r_gone <- run_predict(fx_ab, d_gone, "--min_cov", "80")
check("--min_cov is no longer an option (the run is refused)", r_gone$status != 0L)

# ── --redo yes defeats both caches ────────────────────────────────────────────

d_cache <- new_dir("cache")
invisible(run_predict(fx_ab, d_cache, "--iddef", "2"))
r_stale <- run_predict(fx_ab, d_cache, "--iddef", "1")
check("no --redo: cached F is reused, whatever --iddef made it (existing behaviour)",
      f_at(d_cache) == 0.6667)
check("no --redo: log warns that cached F-measures are reused",
      any(grepl("without --redo yes, F-measures cached", r_stale$log, fixed = TRUE)))
invisible(run_predict(fx_ab, d_cache, "--iddef", "1", "--redo", "yes"))
check("--redo yes: F recomputed under --iddef 1, F = 1", f_at(d_cache) == 1)

# Two global calls with mostly different sequences but the same --out and
# --prefix, as 07a/07c do across ranks. d_def's t.sim stands in for the first
# call's output. Only a_1 is shared with it, so reusing it leaves the two new
# species as singletons: (4 x 0.4 + 4 x 0.4 + 1 x 1) / 9 = 0.4667.
fx_next <- write_fixture(c(A[1], species_group(random_seq(200), "c"),
                           species_group(random_seq(200), "d")),
                         c("Testgenus_alpha", rep(c("Testgenus_gamma", "Testgenus_delta"), each = 4)),
                         file.path(root, "next"))

d_nore <- new_dir("reuse_no_redo")
d_redo <- new_dir("reuse_redo")
invisible(file.copy(file.path(d_def, "t.sim"), d_nore))
invisible(file.copy(file.path(d_def, "t.sim"), d_redo))

invisible(run_predict(fx_next, d_nore, "--iddef", "1"))
check("no --redo: a leftover .sim is reused, F = 0.4667 (existing behaviour)",
      f_at(d_nore) == 0.4667)

r_redo <- run_predict(fx_next, d_redo, "--iddef", "1", "--redo", "yes")
check("--redo yes: leftover .sim is not loaded",
      !any(grepl("Loading similarity matrix", r_redo$log)))
check("--redo yes: similarity recomputed for these sequences, F = 1", f_at(d_redo) == 1)

# ── compute_sim.R, and an explicit --sim under --redo yes ─────────────────────

run_sim <- function(out, ...) {
  system2("Rscript", c("R/compute_sim.R", "--input", fx_ab$fasta, "--out", out,
                       "--n_cpus", "1", "--tmp_dir", file.path(out, "tmp"), ...),
          stdout = FALSE, stderr = FALSE)
}
d_sim <- new_dir("sim")
check("compute_sim.R --iddef 1 runs", run_sim(d_sim, "--iddef", "1") == 0L)
sim <- fread(file.path(d_sim, "ab.sim"), header = FALSE, col.names = c("i", "j", "score"))
ab  <- sim[substr(i, 1, 1) != substr(j, 1, 1)]
check("compute_sim.R --iddef 1: every A-B score is under 0.2", nrow(ab) == 32L && all(ab$score < 0.2))
check("compute_sim.R: all 8 sequences present", setequal(unique(sim$i), names(c(A, B))))

d_sim2 <- new_dir("sim2")
check("compute_sim.R --iddef 2 runs", run_sim(d_sim2, "--iddef", "2") == 0L)
sim2 <- fread(file.path(d_sim2, "ab.sim"), header = FALSE, col.names = c("i", "j", "score"))
check("compute_sim.R --iddef 2: the A-B artefact scores 1.0",
      all(sim2[substr(i, 1, 1) != substr(j, 1, 1), score] == 1))
check("compute_sim.R --min_cov is no longer an option", run_sim(new_dir("sim3"), "--min_cov", "80") != 0L)

d_explicit <- new_dir("explicit")
r_explicit <- run_predict(fx_ab, d_explicit, "--iddef", "1", "--sim", file.path(d_sim, "ab.sim"),
                          "--redo", "yes")
check("--redo yes: an explicit --sim is still loaded",
      any(grepl("Loading similarity matrix", r_explicit$log)))
check("a loaded .sim: log notes that --iddef does not apply to it",
      any(grepl("--iddef and --min_overlap do not apply to a loaded matrix", r_explicit$log, fixed = TRUE)))
check(".sim from compute_sim.R gives the same result, F = 1", f_at(d_explicit) == 1)

# ── Validation ────────────────────────────────────────────────────────────────

d_bad <- new_dir("bad")
for (bad in c("5", "-1")) {
  r_bad <- run_predict(fx_ab, d_bad, paste0("--iddef=", bad))
  check(sprintf("--iddef %s is rejected", bad), r_bad$status != 0L)
}

# ── Result ────────────────────────────────────────────────────────────────────

unlink(root, recursive = TRUE)
cat(sprintf("\nRESULT: %d passed, %d failed\n", n_pass, n_fail))
if (n_fail > 0) quit(status = 1)
