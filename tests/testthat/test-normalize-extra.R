# Additional tests for R/normalize.R covering the validation branches and the
# edgeR (TMM) / DESeq2 (RLE) paths when those packages are available.

test_that("TPM errors when gene_lengths length mismatches rows", {
  counts <- matrix(1:4, nrow = 2,
                   dimnames = list(c("g1", "g2"), c("s1", "s2")))
  expect_error(
    bb_normalize(counts, method = "tpm", gene_lengths = c(1000, 2000, 3000)),
    "must match number of rows"
  )
})

test_that("bb_normalize rejects an unknown method", {
  counts <- matrix(1:4, nrow = 2,
                   dimnames = list(c("g1", "g2"), c("s1", "s2")))
  expect_error(bb_normalize(counts, method = "quantile"), "should be one of")
})

test_that("bb_normalize accepts a data.frame and coerces to matrix", {
  df <- data.frame(s1 = c(10, 20), s2 = c(30, 40),
                   row.names = c("g1", "g2"))
  cpm <- bb_normalize(df, method = "cpm")
  expect_true(is.matrix(cpm))
  expect_equal(unname(colSums(cpm)), c(1e6, 1e6))
})

test_that("TMM normalization works when edgeR is available", {
  skip_if_not_installed("edgeR")
  set.seed(1)
  counts <- matrix(rpois(40, 100), nrow = 10,
                   dimnames = list(paste0("g", 1:10), paste0("s", 1:4)))
  tmm <- bb_normalize(counts, method = "tmm")
  expect_equal(dim(tmm), dim(counts))
  expect_true(all(is.finite(tmm)))
})

test_that("RLE normalization works when DESeq2 is available", {
  skip_if_not_installed("DESeq2")
  set.seed(2)
  counts <- matrix(rpois(60, 100), nrow = 10,
                   dimnames = list(paste0("g", 1:10), paste0("s", 1:6)))
  rle <- bb_normalize(counts, method = "rle")
  expect_equal(dim(rle), dim(counts))
  expect_true(all(rle >= 0))
})
