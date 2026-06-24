# Tests for R/alignment.R
# Focus: input validation, aligner dispatch via match.arg, STAR log parsing,
# and tool-availability checks. The aligners shell out to external binaries
# that are not present in the test environment, so we exercise the pure-R
# validation and parsing branches plus the check_tool() guard.

test_that("bb_align validates paired-end input length", {
  fq <- withr::local_tempfile(fileext = ".fastq")
  writeLines(c("@r1", "ACGT", "+", "IIII"), fq)
  expect_error(
    bb_align(fq, "index", tempfile("aln"), aligner = "STAR", paired = TRUE),
    "length 2"
  )
})

test_that("bb_align errors on missing FASTQ file", {
  expect_error(
    bb_align("/nonexistent/reads.fastq", "index", tempfile("aln"),
             aligner = "STAR"),
    "FASTQ file not found"
  )
})

test_that("bb_align creates output_dir and dispatches to the chosen aligner", {
  fq <- withr::local_tempfile(fileext = ".fastq")
  writeLines(c("@r1", "ACGT", "+", "IIII"), fq)
  out_dir <- file.path(tempdir(), paste0("bb_aln_", as.integer(runif(1, 1, 1e8))))
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  # STAR is not installed: dispatch should reach check_tool() and stop there,
  # which proves switch() routed to .align_star and output_dir was created.
  expect_error(
    bb_align(fq, "index", out_dir, aligner = "STAR"),
    "STAR"
  )
  expect_true(dir.exists(out_dir))
})

test_that("bb_align rejects an unknown aligner", {
  fq <- withr::local_tempfile(fileext = ".fastq")
  writeLines(c("@r1", "ACGT", "+", "IIII"), fq)
  expect_error(
    bb_align(fq, "index", tempfile("aln"), aligner = "bowtie"),
    "should be one of"
  )
})

test_that(".parse_star_log returns status row when log missing", {
  res <- bambamR:::.parse_star_log("/nonexistent/Log.final.out")
  expect_s3_class(res, "data.frame")
  expect_equal(res$metric, "status")
  expect_equal(res$value, "log not found")
})

test_that(".parse_star_log extracts known metrics", {
  log_file <- withr::local_tempfile(fileext = ".out")
  writeLines(c(
    "                          Number of input reads |\t1000",
    "                   Uniquely mapped reads number |\t950",
    "                        Uniquely mapped reads % |\t95.00%",
    "                       Number of splices: Total |\t120",
    "                    Mismatch rate per base, % |\t0.20%"
  ), log_file)
  res <- bambamR:::.parse_star_log(log_file)
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 5L)
  expect_true("Number of input reads" %in% res$metric)
  expect_equal(res$value[res$metric == "Number of input reads"], "1000")
})

test_that(".parse_star_log returns empty frame when no metrics match", {
  log_file <- withr::local_tempfile(fileext = ".out")
  writeLines(c("unrelated line", "another line"), log_file)
  res <- bambamR:::.parse_star_log(log_file)
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 0L)
})
