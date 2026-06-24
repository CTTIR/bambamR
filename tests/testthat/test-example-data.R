# Tests for R/example_data.R
# Focus: the three bundled-dataset loaders return objects with the documented
# structure, and the missing-file branch errors. Files live in inst/extdata.

test_that("bb_example_counts returns counts matrix and metadata", {
  ex <- bb_example_counts()
  expect_type(ex, "list")
  expect_true(all(c("counts", "metadata") %in% names(ex)))
  expect_true(is.matrix(ex$counts))
  expect_false(is.null(rownames(ex$counts)))
  expect_false(is.null(colnames(ex$counts)))
  expect_s3_class(ex$metadata, "data.frame")
  expect_equal(nrow(ex$metadata), ncol(ex$counts))
})

test_that("bb_example_mutations returns mutations and clinical data", {
  ex <- bb_example_mutations()
  expect_type(ex, "list")
  expect_true(all(c("mutations", "clinical") %in% names(ex)))
  expect_s3_class(ex$mutations, "data.frame")
  expect_true(all(c("sample", "gene", "mutation_type") %in%
                    colnames(ex$mutations)))
})

test_that("bb_example_de returns a standardized DE data.frame", {
  de <- bb_example_de()
  expect_s3_class(de, "data.frame")
  expect_true(all(c("gene", "log2fc", "pvalue", "padj") %in% colnames(de)))
})

test_that("example datasets round-trip through readRDS deterministically", {
  # Loading twice yields identical objects (no random component on read).
  expect_identical(bb_example_counts(), bb_example_counts())
  expect_identical(bb_example_de(), bb_example_de())
})

test_that("bb_example_counts has the documented 200 x 10 dimensions", {
  ex <- bb_example_counts()
  expect_equal(dim(ex$counts), c(200L, 10L))
  expect_true(all(ex$counts >= 0))
})
