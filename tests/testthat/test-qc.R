test_that("bb_qc works with FASTQ input", {
  tmp <- tempfile(fileext = ".fastq")
  on.exit(unlink(tmp))

  writeLines(c(
    "@r1", "ACGTACGT", "+", "IIIIIIII",
    "@r2", "GCGCGCGC", "+", "HHHHHHHH",
    "@r3", "AATTCCGG", "+", "FFFFFFFF"
  ), tmp)

  qc <- bb_qc(fastq_path = tmp)
  expect_s3_class(qc, "bb_qc")
  expect_equal(qc$read_counts[[basename(tmp)]], 3)
  expect_true(length(qc$gc_content[[basename(tmp)]]) == 3)
})

test_that("bb_qc_summary returns data.frame", {
  tmp <- tempfile(fileext = ".fastq")
  on.exit(unlink(tmp))

  writeLines(c("@r1", "ACGT", "+", "IIII"), tmp)
  qc <- bb_qc(fastq_path = tmp)
  summ <- bb_qc_summary(qc)
  expect_s3_class(summ, "data.frame")
  expect_true("total_reads" %in% colnames(summ))
})

test_that("bb_qc errors without input", {
  expect_error(bb_qc(), "at least one")
})

test_that("print.bb_qc works", {
  tmp <- tempfile(fileext = ".fastq")
  on.exit(unlink(tmp))
  writeLines(c("@r1", "ACGT", "+", "IIII"), tmp)
  qc <- bb_qc(fastq_path = tmp)
  expect_output(print(qc), "bambamR QC Summary")
})
