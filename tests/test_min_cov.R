#!/usr/bin/env Rscript
# test_min_cov.R — Regression test for --min_cov and --redo yes in R/predict.R
# and R/compute_sim.R.
#
# --iddef 2 excludes terminal gaps, so two unrelated sequences can report 100%
# identity from a few overlapping bases at their ends. The fixture builds that
# case on purpose: every sequence of species A ends with the 30 bases that every
# sequence of species B starts with. Within a species, sequences differ by one
# substitution from a common parent (at least 99% identical to each other).
#
# So at threshold 0.99 the correct clustering is two clusters and F = 1, while a
# spurious A-B link merges them and gives F = 2 x 4 / (4 + 8) = 0.6667.
#
# Runs predict.R as a subprocess, because sourcing it would trigger its CLI
# parsing.
#
# Usage: Rscript tests/test_min_cov.R
# Note:  Must be run from the project root directory, with vsearch on PATH.

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

root <- tempfile("test_min_cov_")
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

# ── The fixture must actually trigger the artefact, or the rest proves nothing ─

hits <- file.path(root, "hits.txt")
system2("vsearch", c("--allpairs_global", fx_ab$fasta, "--acceptall", "--userout", hits,
                     "--userfields", "query+target+id+qcov+tcov", "--iddef", "2",
                     "--threads", "1"), stdout = FALSE, stderr = FALSE)
pairs <- fread(hits, header = FALSE, col.names = c("q", "t", "id", "qcov", "tcov"))
cross <- pairs[substr(q, 1, 1) != substr(t, 1, 1)]
check("fixture: every A-B pair is a 100% match on under 20% coverage",
      nrow(cross) == 16L && all(cross$id == 100) && all(cross$qcov < 20) && all(cross$tcov < 20))

# ── Gate: global prediction ───────────────────────────────────────────────────

d_off <- new_dir("off")
r_off <- run_predict(fx_ab, d_off)
check("default run succeeds", r_off$status == 0L)
check("default (no gate): spurious pairs merge the species, F = 0.6667", f_at(d_off) == 0.6667)

d_on <- new_dir("on")
r_on <- run_predict(fx_ab, d_on, "--min_cov", "80")
check("--min_cov 80 run succeeds", r_on$status == 0L)
check("--min_cov 80: spurious pairs removed, F = 1", f_at(d_on) == 1)

# ── Gate: local prediction, which goes through process_dataset() and workers ──

d_local <- new_dir("local")
r_local <- run_predict(fx_ab, d_local, "--min_cov", "80", higher_rank = "genus")
check("local run with --min_cov 80 succeeds", r_local$status == 0L)
check("local: gate reaches per-dataset vsearch, F = 1",
      f_at(d_local, dataset = "Testgenus") == 1)

d_par <- new_dir("parallel")
r_par <- run_predict(fx_ab, d_par, "--min_cov", "80", higher_rank = "genus",
                     "--run_parallel", "yes", "--n_cpus", "2")
check("parallel local run with --min_cov 80 succeeds", r_par$status == 0L)
check("parallel: gate reaches the workers, F = 1",
      f_at(d_par, dataset = "Testgenus") == 1)

# ── Every pair gated out: sequences stay as singletons, the dataset still scores ─

fx_two <- write_fixture(c(A[1], B[1]), c("Testgenus_alpha", "Testgenus_beta"),
                        file.path(root, "two"))
d_two <- new_dir("two")
r_two <- run_predict(fx_two, d_two, "--min_cov", "80")
check("all pairs gated out: run succeeds", r_two$status == 0L)
check("all pairs gated out: two singletons, F = 1", f_at(d_two, threshold = "1.0000") == 1)

# ── --redo yes defeats both caches ────────────────────────────────────────────

d_cache <- new_dir("cache")
invisible(run_predict(fx_ab, d_cache))
r_stale <- run_predict(fx_ab, d_cache, "--min_cov", "80")
check("no --redo: cached F is reused (existing behaviour)", f_at(d_cache) == 0.6667)
check("no --redo: log warns that cached F-measures are reused",
      any(grepl("--min_cov without --redo yes", r_stale$log)))
invisible(run_predict(fx_ab, d_cache, "--min_cov", "80", "--redo", "yes"))
check("--redo yes: F recomputed under the gate, F = 1", f_at(d_cache) == 1)

# Two global calls with mostly different sequences but the same --out and
# --prefix, as 07a/07c do across ranks. d_on's t.sim stands in for the first
# call's output. Only a_1 is shared with it, so reusing it leaves the two new
# species as singletons: (4 x 0.4 + 4 x 0.4 + 1 x 1) / 9 = 0.4667.
fx_next <- write_fixture(c(A[1], species_group(random_seq(200), "c"),
                           species_group(random_seq(200), "d")),
                         c("Testgenus_alpha", rep(c("Testgenus_gamma", "Testgenus_delta"), each = 4)),
                         file.path(root, "next"))

d_nore <- new_dir("reuse_no_redo")
d_redo <- new_dir("reuse_redo")
invisible(file.copy(file.path(d_on, "t.sim"), d_nore))
invisible(file.copy(file.path(d_on, "t.sim"), d_redo))

invisible(run_predict(fx_next, d_nore))
check("no --redo: a leftover .sim is reused, F = 0.4667 (existing behaviour)",
      f_at(d_nore) == 0.4667)

r_redo <- run_predict(fx_next, d_redo, "--min_cov", "80", "--redo", "yes")
check("--redo yes: leftover .sim is not loaded",
      !any(grepl("Loading similarity matrix", r_redo$log)))
check("--redo yes: similarity recomputed for these sequences, F = 1", f_at(d_redo) == 1)

# ── compute_sim.R, and an explicit --sim under --redo yes ─────────────────────

d_sim <- new_dir("sim")
sim_status <- system2("Rscript",
                      c("R/compute_sim.R", "--input", fx_ab$fasta, "--out", d_sim,
                        "--min_cov", "80", "--n_cpus", "1", "--tmp_dir", file.path(d_sim, "tmp")),
                      stdout = FALSE, stderr = FALSE)
check("compute_sim.R runs", sim_status == 0L)
sim <- fread(file.path(d_sim, "ab.sim"), header = FALSE, col.names = c("i", "j", "score"))
check("compute_sim.R --min_cov 80: no A-B pairs", !any(substr(sim$i, 1, 1) != substr(sim$j, 1, 1)))
check("compute_sim.R --min_cov 80: all 8 sequences present",
      setequal(unique(sim$i), names(c(A, B))))

d_explicit <- new_dir("explicit")
r_explicit <- run_predict(fx_ab, d_explicit, "--sim", file.path(d_sim, "ab.sim"),
                          "--min_cov", "80", "--redo", "yes")
check("--redo yes: an explicit --sim is still loaded",
      any(grepl("Loading similarity matrix", r_explicit$log)))
check("--min_cov with a loaded .sim: log warns it is not applied",
      any(grepl("does not apply to the loaded similarity matrix", r_explicit$log)))
check("gated .sim from compute_sim.R gives the same result, F = 1", f_at(d_explicit) == 1)

# ── Validation ────────────────────────────────────────────────────────────────

d_bad <- new_dir("bad")
for (bad in c("0.8", "101", "-1")) {
  r_bad <- run_predict(fx_ab, d_bad, paste0("--min_cov=", bad))
  check(sprintf("--min_cov %s is rejected", bad), r_bad$status != 0L)
}

# ── Result ────────────────────────────────────────────────────────────────────

unlink(root, recursive = TRUE)
cat(sprintf("\nRESULT: %d passed, %d failed\n", n_pass, n_fail))
if (n_fail > 0) quit(status = 1)
