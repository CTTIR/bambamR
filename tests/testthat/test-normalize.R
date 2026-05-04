test_that("CPM normalization works correctly", {
  counts <- matrix(
    c(10, 20, 30, 40, 50, 60),
    nrow = 3, ncol = 2,
    dimnames = list(paste0("g", 1:3), paste0("s", 1:2))
  )
  cpm <- bb_normalize(counts, method = "cpm")
  expect_equal(dim(cpm), dim(counts))
  # CPM columns should sum to 1e6

  expect_equal(colSums(cpm), c(s1 = 1e6, s2 = 1e6))
})

test_that("TPM normalization works with gene lengths", {
  counts <- matrix(
    c(100, 200, 300, 400),
    nrow = 2, ncol = 2,
    dimnames = list(c("g1", "g2"), c("s1", "s2"))
  )
  gene_lengths <- c(1000, 2000)
  tpm <- bb_normalize(counts, method = "tpm", gene_lengths = gene_lengths)
  expect_equal(dim(tpm), dim(counts))
  # TPM columns should sum to 1e6
  expect_equal(colSums(tpm), c(s1 = 1e6, s2 = 1e6), tolerance = 1)
})

test_that("TPM errors without gene lengths", {
  counts <- matrix(1:4, nrow = 2, dimnames = list(c("g1", "g2"), c("s1", "s2")))
  expect_error(bb_normalize(counts, method = "tpm"),
               "gene_lengths.*required")
})

test_that("normalize validates count matrix", {
  expect_error(bb_normalize("not a matrix", method = "cpm"))
  expect_error(bb_normalize(matrix(c("a", "b"), nrow = 1,
                                    dimnames = list("g1", c("s1", "s2"))),
                             method = "cpm"),
               "numeric")
})

test_that("TMM normalization errors without edgeR", {
  skip_if(requireNamespace("edgeR", quietly = TRUE),
          "edgeR is installed, cannot test fallback")
  counts <- matrix(rpois(20, 100), nrow = 5,
                   dimnames = list(paste0("g", 1:5), paste0("s", 1:4)))
  expect_error(bb_normalize(counts, method = "tmm"), "edgeR")
})
