# Tests for R/counting.R
# Focus: input validation, method = match.arg dispatch, sample-name cleaning,
# and the featureCounts tool guard. Bioconductor counting needs a real BAM +
# GTF, so we validate the guards and the deterministic helper.

test_that("bb_count_reads errors on missing BAM files", {
  ann <- withr::local_tempfile(fileext = ".gtf")
  writeLines("placeholder", ann)
  expect_error(
    bb_count_reads(c("/nope/a.bam", "/nope/b.bam"), ann),
    "BAM file\\(s\\) not found"
  )
})

test_that("bb_count_reads errors on missing annotation", {
  bam <- withr::local_tempfile(fileext = ".bam")
  writeLines("x", bam)
  expect_error(
    bb_count_reads(bam, "/nonexistent/genes.gtf"),
    "Annotation file not found"
  )
})

test_that("bb_count_reads featureCounts method requires the external tool", {
  bam <- withr::local_tempfile(fileext = ".bam")
  writeLines("x", bam)
  ann <- withr::local_tempfile(fileext = ".gtf")
  writeLines("x", ann)
  # featureCounts is not on PATH in CI; check_tool() must stop.
  expect_error(
    bb_count_reads(bam, ann, method = "featureCounts"),
    "featureCounts"
  )
})

test_that("bb_count_reads rejects an unknown method", {
  bam <- withr::local_tempfile(fileext = ".bam")
  writeLines("x", bam)
  ann <- withr::local_tempfile(fileext = ".gtf")
  writeLines("x", ann)
  expect_error(
    bb_count_reads(bam, ann, method = "htseq"),
    "should be one of"
  )
})

test_that(".clean_sample_names strips .bam extension case-insensitively", {
  paths <- c("/data/sampleA.bam", "dir/sampleB.BAM", "x/sampleC.Bam")
  expect_equal(
    bambamR:::.clean_sample_names(paths),
    c("sampleA", "sampleB", "sampleC")
  )
})

test_that(".clean_sample_names leaves non-bam names intact aside from basename", {
  expect_equal(bambamR:::.clean_sample_names("/a/b/foo.sorted.bam"),
               "foo.sorted")
})

test_that("bb_count_reads internal method requires GenomicAlignments", {
  skip_if(requireNamespace("GenomicAlignments", quietly = TRUE) &&
            requireNamespace("GenomicRanges", quietly = TRUE),
          "GenomicAlignments present; missing-package guard not reachable")
  bam <- withr::local_tempfile(fileext = ".bam")
  writeLines("x", bam)
  ann <- withr::local_tempfile(fileext = ".gtf")
  writeLines("x", ann)
  expect_error(bb_count_reads(bam, ann, method = "internal"),
               "GenomicAlignments")
})

test_that("bb_count_reads bioc path stops on the GenomicFeatures guard", {
  skip_if_not_installed("Rsamtools")
  skip_if_not_installed("GenomicAlignments")
  skip_if(requireNamespace("GenomicFeatures", quietly = TRUE),
          "GenomicFeatures present; the guard would not trigger")
  dir <- withr::local_tempdir()
  bam <- make_test_bam(dir, n_mapped = 1L, n_unmapped = 0L)
  skip_if(is.null(bam), "could not build a test BAM")
  ann <- file.path(dir, "genes.gtf")
  writeLines(paste("chr1", "src", "exon", "90", "200", ".", "+", ".",
                   "gene_id \"G1\";", sep = "\t"), ann)
  # method = "auto" reaches .count_bioc, which needs GenomicFeatures.
  expect_error(bb_count_reads(bam, ann, method = "auto"), "GenomicFeatures")
})
