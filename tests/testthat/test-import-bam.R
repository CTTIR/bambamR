# Tests for the Rsamtools-backed paths of R/import_bam.R using a real tiny BAM.

test_that("bb_count_bam counts all records via Rsamtools", {
  skip_if_not_installed("Rsamtools")
  dir <- withr::local_tempdir()
  bam <- make_test_bam(dir, n_mapped = 2L, n_unmapped = 1L)
  skip_if(is.null(bam), "could not build a test BAM")
  expect_equal(bb_count_bam(bam), 3L)
})

test_that("bb_read_bam returns a data.frame of alignments via Rsamtools", {
  skip_if_not_installed("Rsamtools")
  dir <- withr::local_tempdir()
  bam <- make_test_bam(dir, n_mapped = 2L, n_unmapped = 1L)
  skip_if(is.null(bam), "could not build a test BAM")
  df <- bb_read_bam(bam, index = TRUE)
  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 3L)
  expect_true(all(c("qname", "flag", "rname", "pos", "mapq", "cigar") %in%
                    colnames(df)))
})

test_that("bb_read_bam re-uses an existing index without re-indexing", {
  skip_if_not_installed("Rsamtools")
  dir <- withr::local_tempdir()
  bam <- make_test_bam(dir, n_mapped = 1L, n_unmapped = 0L)
  skip_if(is.null(bam), "could not build a test BAM")
  Rsamtools::indexBam(bam)  # create .bai up front
  df <- bb_read_bam(bam, index = TRUE)
  expect_equal(nrow(df), 1L)
})
