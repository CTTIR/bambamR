test_that("pipeline from count matrix works (CRAN-only mode)", {
  set.seed(42)
  counts <- matrix(rpois(600, 100), nrow = 100, ncol = 6,
                   dimnames = list(paste0("g", 1:100), paste0("S", 1:6)))
  sample_info <- data.frame(
    condition = factor(rep(c("ctrl", "treat"), each = 3)),
    row.names = paste0("S", 1:6)
  )

  out_dir <- tempfile("bambamR_test_")
  on.exit(unlink(out_dir, recursive = TRUE))

  # Skip DE since it requires Bioconductor
  result <- bb_pipeline(
    count_matrix = counts,
    sample_info = sample_info,
    output_dir = out_dir,
    skip = c("de")
  )
  expect_s3_class(result, "bb_result")
  expect_equal(dim(result$counts), dim(counts))
})

test_that("pipeline validates inputs correctly", {
  expect_error(bb_pipeline(), "Provide one of")

  expect_error(
    bb_pipeline(count_matrix = matrix(1:4, 2, 2,
                                       dimnames = list(c("g1", "g2"),
                                                       c("s1", "s2")))),
    "sample_info.*required"
  )

  expect_error(
    bb_pipeline(fastq_dir = "/nonexistent/dir"),
    "not found"
  )
})

test_that("bb_result print method works", {
  r <- bambamR:::new_bb_result(
    counts = matrix(1:4, 2, 2, dimnames = list(c("g1", "g2"), c("s1", "s2")))
  )
  expect_output(print(r), "bambamR Result")
  expect_output(print(r), "Genes:\\s+2")
})
