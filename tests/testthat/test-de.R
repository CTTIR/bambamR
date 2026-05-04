test_that("bb_deseq2 requires DESeq2 package", {
  skip_if(requireNamespace("DESeq2", quietly = TRUE),
          "DESeq2 is installed")
  counts <- matrix(rpois(60, 100), nrow = 10,
                   dimnames = list(paste0("g", 1:10), paste0("s", 1:6)))
  coldata <- data.frame(
    condition = factor(rep(c("a", "b"), each = 3)),
    row.names = paste0("s", 1:6)
  )
  expect_error(bb_deseq2(counts, coldata), "DESeq2")
})

test_that("bb_edger requires edgeR package", {
  skip_if(requireNamespace("edgeR", quietly = TRUE),
          "edgeR is installed")
  counts <- matrix(rpois(60, 100), nrow = 10,
                   dimnames = list(paste0("g", 1:10), paste0("s", 1:6)))
  group <- rep(c("a", "b"), each = 3)
  expect_error(bb_edger(counts, group), "edgeR")
})

test_that("bb_limma_voom requires limma", {
  skip_if(requireNamespace("limma", quietly = TRUE),
          "limma is installed")
  counts <- matrix(rpois(60, 100), nrow = 10,
                   dimnames = list(paste0("g", 1:10), paste0("s", 1:6)))
  design <- matrix(c(rep(1, 6), rep(0:1, each = 3)), ncol = 2)
  expect_error(bb_limma_voom(counts, design), "limma")
})

test_that("bb_deseq2 returns standardized result when DESeq2 available", {
  skip_if_not_installed("DESeq2")
  set.seed(42)
  counts <- matrix(rpois(600, 100), nrow = 100, ncol = 6,
                   dimnames = list(paste0("g", 1:100), paste0("s", 1:6)))
  coldata <- data.frame(
    condition = factor(rep(c("ctrl", "treat"), each = 3)),
    row.names = paste0("s", 1:6)
  )
  res <- bb_deseq2(counts, coldata)
  expect_s3_class(res, "data.frame")
  expect_true(all(c("gene", "log2fc", "pvalue", "padj") %in% colnames(res)))
})

test_that("bb_edger returns standardized result when edgeR available", {
  skip_if_not_installed("edgeR")
  set.seed(42)
  counts <- matrix(rpois(600, 100), nrow = 100, ncol = 6,
                   dimnames = list(paste0("g", 1:100), paste0("s", 1:6)))
  group <- factor(rep(c("ctrl", "treat"), each = 3))
  res <- bb_edger(counts, group)
  expect_s3_class(res, "data.frame")
  expect_true(all(c("gene", "log2fc", "pvalue", "padj") %in% colnames(res)))
})
