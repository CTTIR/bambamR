test_that("bb_read_fastq reads simple FASTQ", {
  tmp <- tempfile(fileext = ".fastq")
  on.exit(unlink(tmp))

  writeLines(c(
    "@read1", "ACGTACGT", "+", "IIIIIIII",
    "@read2", "TGCATGCA", "+", "HHHHHHHH",
    "@read3", "GGGGCCCC", "+", "FFFFFFFF"
  ), tmp)

  reads <- bb_read_fastq(tmp)
  expect_s3_class(reads, "data.frame")
  expect_equal(nrow(reads), 3)
  expect_equal(colnames(reads), c("id", "sequence", "quality"))
  expect_equal(reads$id[1], "read1")
  expect_equal(reads$sequence[2], "TGCATGCA")
})

test_that("bb_read_fastq respects n argument", {
  tmp <- tempfile(fileext = ".fastq")
  on.exit(unlink(tmp))

  writeLines(c(
    "@r1", "ACGT", "+", "IIII",
    "@r2", "TGCA", "+", "HHHH",
    "@r3", "GGCC", "+", "FFFF"
  ), tmp)

  reads <- bb_read_fastq(tmp, n = 2)
  expect_equal(nrow(reads), 2)
})

test_that("bb_read_fastq errors on missing file", {
  expect_error(bb_read_fastq("/nonexistent/file.fastq"), "not found")
})

test_that("bb_read_bam errors on missing file", {
  expect_error(bb_read_bam("/nonexistent/file.bam"), "not found")
})

test_that("bb_count_bam errors on missing file", {
  expect_error(bb_count_bam("/nonexistent/file.bam"), "not found")
})
