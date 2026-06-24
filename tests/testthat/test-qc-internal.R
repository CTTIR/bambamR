# Tests for the internal QC helpers in R/qc.R and the BAM branch of bb_qc.

test_that(".compute_gc returns the GC fraction per sequence", {
  res <- bambamR:::.compute_gc(c("GCGC", "ATAT", "GGCC"))
  expect_equal(res, c(1, 0, 1))
})

test_that(".compute_qual_summary returns per-position means", {
  # Phred+33: 'I' = 40, 'H' = 39
  summ <- bambamR:::.compute_qual_summary(c("IIII", "HHHH"))
  expect_s3_class(summ, "data.frame")
  expect_equal(nrow(summ), 4L)
  expect_equal(summ$mean_quality, rep(39.5, 4))
})

test_that(".compute_qual_summary handles ragged read lengths with NA padding", {
  summ <- bambamR:::.compute_qual_summary(c("IIII", "II"))
  expect_equal(nrow(summ), 4L)
  # last two positions only seen in the first read (quality 40)
  expect_equal(summ$mean_quality[4], 40)
})

test_that(".compute_qual_summary returns NULL for no input", {
  expect_null(bambamR:::.compute_qual_summary(character(0)))
})

test_that(".compute_mapping_rate falls back to NA without Rsamtools/samtools", {
  skip_if(requireNamespace("Rsamtools", quietly = TRUE),
          "Rsamtools present; fallback path not reachable here")
  skip_if(nchar(Sys.which("samtools")) > 0L, "samtools present")
  expect_true(is.na(bambamR:::.compute_mapping_rate("/x.bam")))
})

test_that("bb_qc warns and skips missing FASTQ files", {
  expect_warning(
    qc <- bb_qc(fastq_path = "/nonexistent/file.fastq"),
    "not found"
  )
  expect_s3_class(qc, "bb_qc")
  expect_equal(length(qc$read_counts), 0L)
})

test_that("bb_qc_summary errors on a non-bb_qc object", {
  expect_error(bb_qc_summary(list(a = 1)), "must be a bb_qc object")
})

test_that("summary.bb_qc dispatches to bb_qc_summary", {
  fq <- withr::local_tempfile(fileext = ".fastq")
  writeLines(c("@r1", "ACGT", "+", "IIII"), fq)
  qc <- bb_qc(fastq_path = fq)
  s <- summary(qc)
  expect_s3_class(s, "data.frame")
  expect_true("median_gc" %in% colnames(s))
})
