# Tests for the base-R / system fallbacks in R/import_fastq.R and R/import_bam.R.
# The exported wrappers prefer Bioconductor (ShortRead / Rsamtools) when
# installed, so to exercise the pure-R fallback parsers deterministically we
# call the internal helpers directly.

test_that(".read_fastq_base parses a plain FASTQ and strips the @ from ids", {
  fq <- withr::local_tempfile(fileext = ".fastq")
  writeLines(c(
    "@read1 some description", "ACGTACGT", "+", "IIIIIIII",
    "@read2", "TGCATGCA", "+", "HHHHHHHH"
  ), fq)
  res <- bambamR:::.read_fastq_base(fq, n = NULL)
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 2L)
  expect_equal(colnames(res), c("id", "sequence", "quality"))
  expect_equal(res$id[1], "read1 some description")
  expect_equal(res$sequence[2], "TGCATGCA")
})

test_that(".read_fastq_base honours the n argument", {
  fq <- withr::local_tempfile(fileext = ".fastq")
  writeLines(c(
    "@r1", "ACGT", "+", "IIII",
    "@r2", "TGCA", "+", "HHHH",
    "@r3", "GGCC", "+", "FFFF"
  ), fq)
  res <- bambamR:::.read_fastq_base(fq, n = 2L)
  expect_equal(nrow(res), 2L)
})

test_that(".read_fastq_base returns an empty frame for an empty file", {
  fq <- withr::local_tempfile(fileext = ".fastq")
  file.create(fq)
  res <- bambamR:::.read_fastq_base(fq, n = NULL)
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 0L)
})

test_that(".read_fastq_base transparently reads gzipped FASTQ", {
  gz <- withr::local_tempfile(fileext = ".fastq.gz")
  con <- gzfile(gz, "w")
  writeLines(c("@r1", "ACGTACGT", "+", "IIIIIIII"), con)
  close(con)
  res <- bambamR:::.read_fastq_base(gz, n = NULL)
  expect_equal(nrow(res), 1L)
  expect_equal(res$sequence[1], "ACGTACGT")
})

test_that(".read_bam_system guards on a missing samtools tool", {
  skip_if(nchar(Sys.which("samtools")) > 0L,
          "samtools present; cannot test the missing-tool guard")
  expect_error(bambamR:::.read_bam_system("/some/file.bam"), "samtools")
})

test_that(".read_bam_rsamtools returns empty frame for a BAM with no records", {
  skip_if_not_installed("Rsamtools")
  # Build a tiny header-only (zero alignment) BAM via Rsamtools' bundled tools.
  hdr <- c("@HD\tVN:1.6\tSO:coordinate", "@SQ\tSN:chr1\tLN:1000")
  sam <- withr::local_tempfile(fileext = ".sam")
  writeLines(hdr, sam)
  bam <- tryCatch(
    Rsamtools::asBam(sam, overwrite = TRUE,
                     destination = sub("\\.sam$", "", sam)),
    error = function(e) NULL
  )
  skip_if(is.null(bam), "could not build a test BAM")
  withr::defer(unlink(c(bam, paste0(bam, ".bai"))))
  res <- bambamR:::.read_bam_rsamtools(bam, index = TRUE,
                                       what = c("qname", "flag", "rname",
                                                "pos", "mapq", "cigar"))
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 0L)
})
