# Tests for the internal helpers in R/utils.R.

test_that("check_pkg returns invisibly TRUE for an installed package", {
  expect_true(bambamR:::check_pkg("stats"))
})

test_that("check_pkg errors with a Bioconductor install hint", {
  expect_error(
    bambamR:::check_pkg("definitelyNotARealPkg123",
                        reason = "for testing"),
    "install.packages"
  )
})

test_that("check_pkg suggests BiocManager for known Bioc packages", {
  # Force the not-installed path by temporarily mocking requireNamespace so the
  # branch builds a BiocManager hint regardless of what is installed.
  local_mocked_bindings(
    requireNamespace = function(pkg, ...) FALSE,
    .package = "base"
  )
  expect_error(bambamR:::check_pkg("DESeq2"), "BiocManager::install")
})

test_that("check_tool errors when a tool is missing from PATH", {
  expect_error(bambamR:::check_tool("totally_missing_tool_xyz"),
               "not found on PATH")
})

test_that("validate_counts enforces matrix-like numeric input with names", {
  expect_error(bambamR:::validate_counts("nope"), "matrix or data.frame")
  expect_error(
    bambamR:::validate_counts(matrix("a", 1, 1,
                                     dimnames = list("g", "s"))),
    "numeric"
  )
  expect_error(
    bambamR:::validate_counts(matrix(1:4, 2, 2,
                                     dimnames = list(NULL, c("s1", "s2")))),
    "rownames"
  )
  expect_error(
    bambamR:::validate_counts(matrix(1:4, 2, 2,
                                     dimnames = list(c("g1", "g2"), NULL))),
    "colnames"
  )
})

test_that("validate_counts coerces a data.frame to a numeric matrix", {
  df <- data.frame(s1 = c(1, 2), s2 = c(3, 4), row.names = c("g1", "g2"))
  m <- bambamR:::validate_counts(df)
  expect_true(is.matrix(m))
  expect_equal(rownames(m), c("g1", "g2"))
})

test_that("CPM normalization output is reproducible (snapshot value)", {
  counts <- matrix(
    c(10, 20, 30, 40, 50, 60),
    nrow = 3, ncol = 2,
    dimnames = list(paste0("g", 1:3), paste0("s", 1:2))
  )
  cpm <- bb_normalize(counts, method = "cpm")
  expect_snapshot_value(round(cpm, 4), style = "json2")
})
