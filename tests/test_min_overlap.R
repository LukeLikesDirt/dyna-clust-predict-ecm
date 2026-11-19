#!/usr/bin/env Rscript
# test_min_overlap.R — Regression test for --iddef 2 with a minimum overlap in
# bp (--min_overlap), in R/predict.R, R/compute_sim.R and R/remove_complexes.R.
# Synced to the sibling repos (tools/sync_manifest.tsv); the per-region values
# in this repo's scripts are checked by tests/test_script_params.R.
#
# Fixture 1 (the artefact the overlap exists for): every sequence of species A
# ends with the 30 bases that every sequence of species B starts with. --iddef 2
# excludes terminal gaps, so it scores each A-B pair 100% from those 30 bases.
# Within a species, sequences differ by one substitution from a common parent
# (at least 99% identical to each other). At threshold 0.99 the correct
# clustering is two clusters and F = 1, while a spurious A-B link merges them
# and gives F = 2 x 4 / (4 + 8) = 0.6667.
#
# Fixture 2 (what the overlap must keep): species Q's sequences are a 150 bp
# fragment of species P's 300 bp sequences, identical over that overlap. They
# score 100% under --iddef 2 on 150 bp, but cover only 50% of the P sequence, so
# a fractional coverage gate (the former --min_cov 80) would have dropped them.
#
# Runs the R scripts as subprocesses, because sourcing them would trigger their
# CLI parsing.
#
# Usage: Rscript tests/test_min_overlap.R
# Note:  Must be run from the project root directory, with vsearch on PATH.

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

root <- tempfile("test_min_overlap_")
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

# ── Fixtures ──────────────────────────────────────────────────────────────────

set.seed(42)
bases      <- c("A", "C", "G", "T")
random_seq <- function(n) paste(sample(bases, n, replace = TRUE), collapse = "")

# Parent plus three variants, each with one substitution at its own position.
species_group <- function(parent, name, positions = c(60L, 100L, 140L)) {
  chars    <- strsplit(parent, "")[[1]]
  variants <- vapply(positions, function(pos) {
    v <- chars
    v[pos] <- setdiff(bases, v[pos])[1]
    paste(v, collapse = "")
  }, character(1))
  setNames(c(parent, variants), paste0(name, "_", 1:4))
}

shared <- random_seq(30)
A <- species_group(paste0(random_seq(170), shared), "a")
B <- species_group(paste0(shared, random_seq(170)), "b")

big   <- random_seq(300)
P     <- species_group(big, "p", positions = c(30L, 150L, 270L))
Q     <- species_group(substr(big, 101, 250), "q", positions = c(20L, 70L, 120L))

write_fixture <- function(seqs, species, stem, genus = "Testgenus") {
  fasta <- paste0(stem, ".fasta")
  cls   <- paste0(stem, ".classification")
  writeLines(paste0(">", names(seqs), "\n", seqs), fasta)
  fwrite(data.table(id = names(seqs), genus = genus, species = species), cls, sep = "\t")
  list(fasta = fasta, cls = cls)
}

fx_ab <- write_fixture(c(A, B), rep(c("Testgenus_alpha", "Testgenus_beta"), each = 4),
                       file.path(root, "ab"))
fx_pq <- write_fixture(c(P, Q), rep(c("Testgenus_pee", "Testgenus_que"), each = 4),
                       file.path(root, "pq"))

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

run_sim <- function(fx, out, ...) {
  status <- system2("Rscript", c("R/compute_sim.R", "--input", fx$fasta, "--out", out,
                                 "--n_cpus", "1", "--tmp_dir", file.path(out, "tmp"), ...),
                    stdout = FALSE, stderr = FALSE)
  status
}
read_sim <- function(out, stem) {
  fread(file.path(out, paste0(stem, ".sim")), header = FALSE, col.names = c("i", "j", "score"))
}

all_pairs <- function(fasta, iddef = 2) {
  hits <- file.path(root, "hits.txt")
  system2("vsearch", c("--allpairs_global", fasta, "--acceptall", "--userout", hits,
                       "--userfields", "query+target+id+alnlen+qcov+tcov", "--iddef", iddef,
                       "--threads", "1"), stdout = FALSE, stderr = FALSE)
  fread(hits, header = FALSE, col.names = c("q", "t", "id", "alnlen", "qcov", "tcov"))
}

# ── The fixtures must actually trigger what they test ─────────────────────────

cross <- all_pairs(fx_ab$fasta)[substr(q, 1, 1) != substr(t, 1, 1)]
check("fixture 1: every A-B pair is a 100% match over exactly 30 bp",
      nrow(cross) == 16L && all(cross$id == 100) && all(cross$alnlen == 30))
pq <- all_pairs(fx_pq$fasta)
pq <- pq[substr(q, 1, 1) != substr(t, 1, 1)]
check("fixture 2: the Q fragment matches P at 100% over 150 bp on 50% of the P sequence",
      any(pq$id == 100 & pq$alnlen == 150 & (pq$qcov == 50 | pq$tcov == 50)))

# ── predict.R, global ─────────────────────────────────────────────────────────

d_def <- new_dir("default")
r_def <- run_predict(fx_ab, d_def)
check("default run succeeds", r_def$status == 0L)
check("default is --iddef 2 with no overlap filter: the artefact merges A and B, F = 0.6667",
      f_at(d_def) == 0.6667)
check("log notes that --iddef 2 without --min_overlap lets a few shared bases score 100%",
      any(grepl("--min_overlap is 0", r_def$log, fixed = TRUE)))

d_50 <- new_dir("overlap50")
r_50 <- run_predict(fx_ab, d_50, "--min_overlap", "50")
check("--min_overlap 50 run succeeds", r_50$status == 0L)
check("--min_overlap 50: the 30 bp artefact is dropped, F = 1", f_at(d_50) == 1)
check("--min_overlap 50: no --min_overlap note in the log",
      !any(grepl("--min_overlap is 0", r_50$log, fixed = TRUE)))

d_30 <- new_dir("overlap30")
invisible(run_predict(fx_ab, d_30, "--min_overlap", "30"))
check("--min_overlap 30: a pair overlapping by exactly 30 bp is kept (>=), F = 0.6667",
      f_at(d_30) == 0.6667)
d_31 <- new_dir("overlap31")
invisible(run_predict(fx_ab, d_31, "--min_overlap", "31"))
check("--min_overlap 31: a pair overlapping by 30 bp is dropped, F = 1", f_at(d_31) == 1)

# A fragment identical to part of a longer sequence is kept by a bp threshold at
# or below its length, however small a share of the long sequence it covers.
# At 0.99 P and Q are then one cluster of 8, so F = 2 x 4 / (4 + 8) = 0.6667. (A
# threshold above 150 would also reject the Q sequences' own 150 bp overlaps, so
# the dropped case is checked on the similarity matrix, below, not on F.)
d_pq100 <- new_dir("pq100")
invisible(run_predict(fx_pq, d_pq100, "--min_overlap", "100"))
check("a 150 bp fragment of a 300 bp sequence (50% coverage) is kept at --min_overlap 100, F = 0.6667",
      f_at(d_pq100) == 0.6667)

# ── predict.R, local and parallel (the workers compute their own similarity) ──

d_local <- new_dir("local")
r_local <- run_predict(fx_ab, d_local, "--min_overlap", "50", higher_rank = "genus")
check("local run succeeds", r_local$status == 0L)
check("local: per-dataset vsearch honours --min_overlap, F = 1",
      f_at(d_local, dataset = "Testgenus") == 1)

d_local0 <- new_dir("local0")
invisible(run_predict(fx_ab, d_local0, higher_rank = "genus"))
check("local without --min_overlap: the artefact merges A and B, F = 0.6667",
      f_at(d_local0, dataset = "Testgenus") == 0.6667)

d_par <- new_dir("parallel")
r_par <- run_predict(fx_ab, d_par, "--min_overlap", "50", higher_rank = "genus",
                     "--run_parallel", "yes", "--n_cpus", "2")
check("parallel local run succeeds", r_par$status == 0L)
check("parallel: --min_overlap reaches the workers, F = 1",
      f_at(d_par, dataset = "Testgenus") == 1)

# ── A sequence with every pair rejected must stay in its dataset ──────────────
# At --min_overlap 250 every pair is rejected (the longest overlap is 200 bp),
# so vsearch writes an empty file and the 8 sequences are self-only nodes. Each
# is its own cluster, so each species of 4 has best Dice 2 x 1 / (4 + 1) = 0.4
# and F = 0.4. A sequence that dropped out of the matrix would shrink its class.

d_all <- new_dir("allrejected")
r_all <- run_predict(fx_ab, d_all, "--min_overlap", "250")
check("every pair rejected: run succeeds (empty vsearch output is not an error)",
      r_all$status == 0L)
check("every pair rejected: all 8 sequences stay as singleton clusters, F = 0.4",
      f_at(d_all) == 0.4)

# ── compute_sim.R ─────────────────────────────────────────────────────────────

d_sim0 <- new_dir("sim0")
check("compute_sim.R runs (default --iddef 2, no overlap)", run_sim(fx_ab, d_sim0) == 0L)
sim0 <- read_sim(d_sim0, "ab")
check("compute_sim.R default: the A-B artefact scores 1.0",
      all(sim0[substr(i, 1, 1) != substr(j, 1, 1), score] == 1))

d_sim50 <- new_dir("sim50")
check("compute_sim.R --min_overlap 50 runs", run_sim(fx_ab, d_sim50, "--min_overlap", "50") == 0L)
sim50 <- read_sim(d_sim50, "ab")
check("compute_sim.R --min_overlap 50: no A-B pair is written",
      nrow(sim50[substr(i, 1, 1) != substr(j, 1, 1)]) == 0L)
check("compute_sim.R --min_overlap 50: within-species pairs are unchanged",
      nrow(sim50[substr(i, 1, 1) == substr(j, 1, 1) & i != j]) ==
        nrow(sim0[substr(i, 1, 1) == substr(j, 1, 1) & i != j]))
check("compute_sim.R --min_overlap 50: all 8 sequences present", setequal(unique(sim50$i), names(c(A, B))))

d_simall <- new_dir("simall")
check("compute_sim.R with every pair rejected runs", run_sim(fx_ab, d_simall, "--min_overlap", "250") == 0L)
simall <- read_sim(d_simall, "ab")
check("compute_sim.R with every pair rejected: 8 self-only rows",
      nrow(simall) == 8L && all(simall$i == simall$j) && all(simall$score == 1))

d_simpq <- new_dir("simpq")
invisible(run_sim(fx_pq, d_simpq, "--min_overlap", "100"))
simpq <- read_sim(d_simpq, "pq")
check("compute_sim.R keeps the 150 bp fragment link at --min_overlap 100",
      nrow(simpq[substr(i, 1, 1) != substr(j, 1, 1) & score == 1]) > 0L)
d_simpq2 <- new_dir("simpq2")
invisible(run_sim(fx_pq, d_simpq2, "--min_overlap", "200"))
check("compute_sim.R drops it at --min_overlap 200",
      nrow(read_sim(d_simpq2, "pq")[substr(i, 1, 1) != substr(j, 1, 1) & score == 1]) == 0L)

# ── An explicit --sim from compute_sim.R is loaded, with a note ───────────────

d_explicit <- new_dir("explicit")
r_explicit <- run_predict(fx_ab, d_explicit, "--sim", file.path(d_sim50, "ab.sim"), "--redo", "yes")
check("an explicit --sim is loaded", any(grepl("Loading similarity matrix", r_explicit$log)))
check("a loaded .sim: log notes that --iddef and --min_overlap do not apply to it",
      any(grepl("--iddef and --min_overlap do not apply to a loaded matrix", r_explicit$log, fixed = TRUE)))
check(".sim from compute_sim.R --min_overlap 50 gives the same result, F = 1", f_at(d_explicit) == 1)

# ── remove_complexes.R ────────────────────────────────────────────────────────

run_rc <- function(fx, out, ...) {
  dir.create(out, showWarnings = FALSE)
  args <- c("R/remove_complexes.R", "--fasta_in", fx$fasta, "--classification_in", fx$cls,
            "--fasta_out", file.path(out, "out.fasta"),
            "--classification_out", file.path(out, "out.tsv"),
            "--manifest_out", file.path(out, "manifest.txt"),
            "--threshold", "1.0", "--require_complete", "no",
            "--run_parallel", "no", "--n_cpus", "1", "--tmp_dir", file.path(out, "tmp"), ...)
  log    <- suppressWarnings(system2("Rscript", args, stdout = TRUE, stderr = TRUE))
  status <- attr(log, "status")
  list(status = if (is.null(status)) 0L else status, log = log)
}
n_kept <- function(out) {
  length(grep("^>", readLines(file.path(out, "out.fasta"))))
}

rc0 <- run_rc(fx_ab, new_dir("rc_ab_default"))
check("remove_complexes.R default run succeeds", rc0$status == 0L)
check("remove_complexes.R default (--iddef 2, no overlap): the 30 bp artefact makes a complex, 4 of 8 kept",
      n_kept(file.path(root, "rc_ab_default")) == 4L)

rc50 <- run_rc(fx_ab, new_dir("rc_ab_50"), "--min_overlap", "50")
check("remove_complexes.R --min_overlap 50: the artefact is not a complex, all 8 kept",
      rc50$status == 0L && n_kept(file.path(root, "rc_ab_50")) == 8L)

rc1 <- run_rc(fx_ab, new_dir("rc_ab_iddef1"), "--iddef", "1")
check("remove_complexes.R --iddef 1: nothing is identical, all 8 kept",
      rc1$status == 0L && n_kept(file.path(root, "rc_ab_iddef1")) == 8L)

rcpq <- run_rc(fx_pq, new_dir("rc_pq_100"), "--min_overlap", "100")
check("remove_complexes.R --min_overlap 100: a species identical to a 150 bp part of another is a complex, 4 of 8 kept",
      rcpq$status == 0L && n_kept(file.path(root, "rc_pq_100")) == 4L)
rcpq2 <- run_rc(fx_pq, new_dir("rc_pq_200"), "--min_overlap", "200")
check("remove_complexes.R --min_overlap 200: the 150 bp overlap is too short, all 8 kept",
      rcpq2$status == 0L && n_kept(file.path(root, "rc_pq_200")) == 8L)

rcall <- run_rc(fx_ab, new_dir("rc_ab_250"), "--min_overlap", "250")
check("remove_complexes.R with every pair rejected succeeds and keeps all 8",
      rcall$status == 0L && n_kept(file.path(root, "rc_ab_250")) == 8L)

# ── Validation ────────────────────────────────────────────────────────────────

check("predict.R rejects --min_overlap -1",
      run_predict(fx_ab, new_dir("bad1"), "--min_overlap=-1")$status != 0L)
check("compute_sim.R rejects --min_overlap -1",
      run_sim(fx_ab, new_dir("bad2"), "--min_overlap=-1") != 0L)
check("remove_complexes.R rejects --min_overlap -1",
      run_rc(fx_ab, new_dir("bad3"), "--min_overlap=-1")$status != 0L)

# ── Result ────────────────────────────────────────────────────────────────────

unlink(root, recursive = TRUE)
cat(sprintf("\nRESULT: %d passed, %d failed\n", n_pass, n_fail))
if (n_fail > 0) quit(status = 1)
