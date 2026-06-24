# Tests for the BAM branch of bb_qc() and the Rsamtools mapping-rate path.

test_that("bb_qc computes read counts and mapping rate from a BAM", {
  skip_if_not_installed("Rsamtools")
  dir <- withr::local_tempdir()
  bam <- make_test_bam(dir, n_mapped = 3L, n_unmapped = 1L)
  skip_if(is.null(bam), "could not build a test BAM")

  qc <- bb_qc(bam_path = bam)
  expect_s3_class(qc, "bb_qc")
  fname <- basename(bam)
  expect_equal(qc$read_counts[[fname]], 4L)
  # 3 of 4 reads are mapped
  expect_equal(qc$mapping_rate[[fname]], 0.75, tolerance = 1e-8)
})

test_that(".compute_mapping_rate returns 0 for an empty BAM", {
  skip_if_not_installed("Rsamtools")
  dir <- withr::local_tempdir()
  bam <- make_test_bam(dir, n_mapped = 0L, n_unmapped = 0L)
  skip_if(is.null(bam), "could not build a test BAM")
  expect_equal(bambamR:::.compute_mapping_rate(bam), 0)
})

test_that("print.bb_qc shows mapping percentage for BAM input", {
  skip_if_not_installed("Rsamtools")
  dir <- withr::local_tempdir()
  bam <- make_test_bam(dir, n_mapped = 2L, n_unmapped = 0L)
  skip_if(is.null(bam), "could not build a test BAM")
  qc <- bb_qc(bam_path = bam)
  expect_output(print(qc), "mapped")
})

test_that("bb_qc warns and skips a missing BAM", {
  expect_warning(qc <- bb_qc(bam_path = "/nonexistent/x.bam"), "not found")
  expect_s3_class(qc, "bb_qc")
})
