# Extra tests for R/de_analysis.R: validation branches and limma-voom /
# DESeq2 contrast paths when the packages are available.

test_that("bb_deseq2 errors with fewer than 2 matching samples", {
  skip_if_not_installed("DESeq2")
  counts <- matrix(rpois(20, 100), nrow = 10,
                   dimnames = list(paste0("g", 1:10), c("s1", "s2")))
  coldata <- data.frame(condition = factor(c("a", "b")),
                        row.names = c("x1", "x2"))  # no overlap
  expect_error(bb_deseq2(counts, coldata), "at least 2 matching samples")
})

test_that("bb_edger errors when group length mismatches columns", {
  skip_if_not_installed("edgeR")
  counts <- matrix(rpois(60, 100), nrow = 10,
                   dimnames = list(paste0("g", 1:10), paste0("s", 1:6)))
  expect_error(bb_edger(counts, group = c("a", "b")),
               "must equal number of columns")
})

test_that("bb_deseq2 honours an explicit contrast", {
  skip_if_not_installed("DESeq2")
  set.seed(11)
  counts <- matrix(rpois(600, 100), nrow = 100, ncol = 6,
                   dimnames = list(paste0("g", 1:100), paste0("s", 1:6)))
  coldata <- data.frame(
    condition = factor(rep(c("ctrl", "treat"), each = 3)),
    row.names = paste0("s", 1:6)
  )
  res <- bb_deseq2(counts, coldata,
                   contrast = c("condition", "treat", "ctrl"))
  expect_s3_class(res, "data.frame")
  expect_true("basemean" %in% colnames(res))
})

test_that("bb_limma_voom returns a standardized result", {
  skip_if_not_installed("limma")
  skip_if_not_installed("edgeR")
  set.seed(12)
  counts <- matrix(rpois(600, 100), nrow = 100, ncol = 6,
                   dimnames = list(paste0("g", 1:100), paste0("s", 1:6)))
  group <- factor(rep(c("ctrl", "treat"), each = 3))
  design <- stats::model.matrix(~ group)
  res <- bb_limma_voom(counts, design)
  expect_s3_class(res, "data.frame")
  expect_true(all(c("gene", "log2fc", "pvalue", "padj") %in% colnames(res)))
})

test_that("validate_de_result rejects frames missing required columns", {
  expect_error(
    bambamR:::validate_de_result(data.frame(gene = "a", log2fc = 1)),
    "missing required columns"
  )
})

test_that("standardize_de_result attaches basemean only when present", {
  df <- data.frame(gene = "a", log2fc = 1, pvalue = 0.1, padj = 0.2,
                   baseMean = 5)
  res <- bambamR:::standardize_de_result(df, basemean_col = "baseMean")
  expect_true("basemean" %in% colnames(res))
  res2 <- bambamR:::standardize_de_result(df, basemean_col = "missing")
  expect_false("basemean" %in% colnames(res2))
})
