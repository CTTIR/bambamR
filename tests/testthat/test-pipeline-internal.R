# Tests for the internal pipeline helpers and extra bb_pipeline branches
# in R/pipeline.R.

test_that(".determine_start picks the earliest provided entry point", {
  expect_equal(bambamR:::.determine_start("fq", NULL, NULL), "fastq")
  expect_equal(bambamR:::.determine_start(NULL, "bam", NULL), "bam")
  expect_equal(bambamR:::.determine_start(NULL, NULL, matrix(1)), "counts")
  expect_error(bambamR:::.determine_start(NULL, NULL, NULL), "Provide one of")
})

test_that(".validate_pipeline_inputs enforces FASTQ requirements", {
  dir <- withr::local_tempdir()
  expect_error(
    bambamR:::.validate_pipeline_inputs("fastq", dir, NULL, NULL,
                                        NULL, "ann.gtf", NULL),
    "genome_index"
  )
  expect_error(
    bambamR:::.validate_pipeline_inputs("fastq", dir, NULL, NULL,
                                        "idx", NULL, NULL),
    "annotation"
  )
  expect_error(
    bambamR:::.validate_pipeline_inputs("fastq", "/no/such/dir", NULL, NULL,
                                        "idx", "ann", NULL),
    "fastq_dir not found"
  )
})

test_that(".validate_pipeline_inputs enforces BAM requirements", {
  dir <- withr::local_tempdir()
  expect_error(
    bambamR:::.validate_pipeline_inputs("bam", NULL, dir, NULL,
                                        NULL, NULL, NULL),
    "annotation"
  )
  expect_error(
    bambamR:::.validate_pipeline_inputs("bam", NULL, "/no/such/dir", NULL,
                                        NULL, "ann", NULL),
    "bam_dir not found"
  )
})

test_that(".list_fastq finds FASTQ files and errors when none exist", {
  dir <- withr::local_tempdir()
  writeLines(c("@r1", "ACGT", "+", "IIII"), file.path(dir, "a.fastq"))
  writeLines(c("@r1", "ACGT", "+", "IIII"), file.path(dir, "b.fq"))
  found <- bambamR:::.list_fastq(dir)
  expect_length(found, 2L)

  empty_dir <- withr::local_tempdir()
  expect_error(bambamR:::.list_fastq(empty_dir), "No FASTQ files found")
})

test_that("bb_pipeline errors when BAM start finds no BAM files", {
  bam_dir <- withr::local_tempdir()  # empty
  ann <- withr::local_tempfile(fileext = ".gtf")
  writeLines("x", ann)
  expect_error(
    bb_pipeline(bam_dir = bam_dir, annotation = ann,
                sample_info = data.frame(condition = "a")),
    "No BAM files found"
  )
})

test_that("bb_pipeline from counts skips DE and still builds plots", {
  set.seed(7)
  counts <- matrix(rpois(600, 100), nrow = 100, ncol = 6,
                   dimnames = list(paste0("g", 1:100), paste0("S", 1:6)))
  sample_info <- data.frame(
    condition = factor(rep(c("ctrl", "treat"), each = 3)),
    row.names = paste0("S", 1:6)
  )
  out_dir <- withr::local_tempdir()
  result <- bb_pipeline(count_matrix = counts, sample_info = sample_info,
                        output_dir = out_dir, skip = c("de"))
  expect_s3_class(result, "bb_result")
  expect_true("pca" %in% names(result$plots))
  expect_true(file.exists(file.path(out_dir, "bambamR_result.rds")))
})

test_that("bb_pipeline can skip visualization entirely", {
  set.seed(8)
  counts <- matrix(rpois(200, 50), nrow = 50, ncol = 4,
                   dimnames = list(paste0("g", 1:50), paste0("S", 1:4)))
  sample_info <- data.frame(
    condition = factor(rep(c("a", "b"), each = 2)),
    row.names = paste0("S", 1:4)
  )
  out_dir <- withr::local_tempdir()
  result <- bb_pipeline(count_matrix = counts, sample_info = sample_info,
                        output_dir = out_dir, skip = c("de", "viz"))
  expect_s3_class(result, "bb_result")
  expect_equal(length(result$plots), 0L)
})

test_that("bb_pipeline runs the DESeq2 DE step and builds volcano/MA plots", {
  skip_if_not_installed("DESeq2")
  set.seed(9)
  counts <- matrix(rpois(600, 100), nrow = 100, ncol = 6,
                   dimnames = list(paste0("g", 1:100), paste0("S", 1:6)))
  sample_info <- data.frame(
    condition = factor(rep(c("ctrl", "treat"), each = 3)),
    row.names = paste0("S", 1:6)
  )
  out_dir <- withr::local_tempdir()
  result <- bb_pipeline(count_matrix = counts, sample_info = sample_info,
                        output_dir = out_dir, de_method = "DESeq2")
  expect_s3_class(result, "bb_result")
  expect_s3_class(result$de_results, "data.frame")
  expect_true("volcano" %in% names(result$plots))
  expect_true("ma" %in% names(result$plots))
})

test_that(".run_de dispatches to edgeR", {
  skip_if_not_installed("edgeR")
  set.seed(10)
  counts <- matrix(rpois(600, 100), nrow = 100, ncol = 6,
                   dimnames = list(paste0("g", 1:100), paste0("S", 1:6)))
  sample_info <- data.frame(
    condition = factor(rep(c("ctrl", "treat"), each = 3)),
    row.names = paste0("S", 1:6)
  )
  res <- bambamR:::.run_de(counts, sample_info, "edgeR", ~ condition)
  expect_s3_class(res, "data.frame")
  expect_true(all(c("gene", "log2fc", "pvalue", "padj") %in% colnames(res)))
})

test_that(".run_de dispatches to limma", {
  skip_if_not_installed("limma")
  skip_if_not_installed("edgeR")
  set.seed(13)
  counts <- matrix(rpois(600, 100), nrow = 100, ncol = 6,
                   dimnames = list(paste0("g", 1:100), paste0("S", 1:6)))
  sample_info <- data.frame(
    condition = factor(rep(c("ctrl", "treat"), each = 3)),
    row.names = paste0("S", 1:6)
  )
  res <- bambamR:::.run_de(counts, sample_info, "limma", ~ condition)
  expect_s3_class(res, "data.frame")
})

test_that(".run_de edgeR path requires a condition column", {
  counts <- matrix(rpois(24, 100), nrow = 6, ncol = 4,
                   dimnames = list(paste0("g", 1:6), paste0("S", 1:4)))
  sample_info <- data.frame(grp = rep(c("a", "b"), each = 2),
                            row.names = paste0("S", 1:4))
  expect_error(
    bambamR:::.run_de(counts, sample_info, "edgeR", ~ grp),
    "must contain a 'condition' column"
  )
})

test_that("print.bb_result reports DE results and plots", {
  r <- bambamR:::new_bb_result(
    counts = matrix(1:4, 2, 2,
                    dimnames = list(c("g1", "g2"), c("s1", "s2"))),
    de_results = data.frame(gene = c("g1", "g2"), log2fc = c(1, -1),
                            pvalue = c(0.01, 0.2), padj = c(0.02, 0.3)),
    plots = list(pca = ggplot2::ggplot())
  )
  expect_output(print(r), "DE results: available")
  expect_output(print(r), "Plots:")
})
