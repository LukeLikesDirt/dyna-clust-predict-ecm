#!/usr/bin/env Rscript
# test_script_params.R — Static checks that this repo's scripts pass the
# similarity settings the shared R/compute_sim.R and R/predict.R expect. Not
# synced from dyna-clust-predict: these scripts are this repo's own.
#
# Usage: Rscript tests/test_script_params.R
# Note:  Must be run from the project root directory.

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

min_len <- sub("^readonly MIN_LENGTH=", "",
               grep("^readonly MIN_LENGTH=", readLines("scripts/02_derep_and_clean.sh"), value = TRUE))

for (sc in c("scripts/05_compute_sim.sh", "scripts/06_predict_ecm_cutoffs.sh")) {
  lines <- readLines(sc)
  name  <- basename(sc)
  check(sprintf("%s: IDDEF=2", name), any(lines == "readonly IDDEF=2"))
  check(sprintf("%s: MIN_OVERLAP equals 02's MIN_LENGTH (%s bp)", name, min_len),
        length(min_len) == 1L && any(lines == paste0("readonly MIN_OVERLAP=", min_len)))
  n_calls <- if (grepl("^05", name)) 1L else 2L
  check(sprintf("%s: --iddef and --min_overlap passed at all %d call site(s)", name, n_calls),
        sum(grepl('--iddef +"\\$IDDEF"', lines)) == n_calls &&
          sum(grepl('--min_overlap +"\\$MIN_OVERLAP"', lines)) == n_calls)
  check(sprintf("%s: does not pass the removed --min_cov", name),
        !any(grepl("min_cov|MIN_COV", lines)))
}

check("06 names predict.R's output <prefix>_raw_cutoffs",
      any(grepl("_raw_cutoffs", readLines("scripts/06_predict_ecm_cutoffs.sh"), fixed = TRUE)))
check("predict.R writes <prefix>_raw_cutoffs.json (and .txt)",
      any(grepl('"_raw_cutoffs.json"', readLines("R/predict.R"), fixed = TRUE)))

cat(sprintf("\nRESULT: %d passed, %d failed\n", n_pass, n_fail))
if (n_fail > 0) quit(status = 1)
