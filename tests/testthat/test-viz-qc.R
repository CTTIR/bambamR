# Tests for R/viz_qc.R: the bb_plot_qc panel builder, each per-metric plot,
# the empty-data placeholders, and the non-patchwork fallback.

make_qc <- function() {
  fq <- withr::local_tempfile(fileext = ".fastq", .local_envir = parent.frame())
  writeLines(c(
    "@r1", "ACGTACGT", "+", "IIIIIIII",
    "@r2", "GCGCGCGC", "+", "HHHHHHHH",
    "@r3", "ATATGCGC", "+", "GGGGFFFF"
  ), fq)
  bb_qc(fastq_path = fq)
}

test_that("bb_plot_qc validates the object class", {
  expect_error(bb_plot_qc(list(x = 1)), "must be a bb_qc object")
})

test_that("bb_plot_qc rejects an unknown 'which' value", {
  qc <- make_qc()
  expect_error(bb_plot_qc(qc, which = "entropy"), "should be one of")
})

test_that("bb_plot_qc returns a single ggplot for each metric", {
  qc <- make_qc()
  expect_s3_class(bb_plot_qc(qc, which = "quality"), "gg")
  expect_s3_class(bb_plot_qc(qc, which = "gc"), "gg")
  expect_s3_class(bb_plot_qc(qc, which = "length"), "gg")
})

test_that("bb_plot_qc which='all' composes a multi-panel object", {
  qc <- make_qc()
  p <- bb_plot_qc(qc, which = "all")
  expect_true(inherits(p, "patchwork") || inherits(p, "gg"))
})

test_that("bb_plot_qc falls back to single plot when patchwork is absent", {
  qc <- make_qc()
  local_mocked_bindings(
    requireNamespace = function(pkg, ...) {
      if (identical(pkg, "patchwork")) FALSE else TRUE
    },
    .package = "base"
  )
  expect_message(p <- bb_plot_qc(qc, which = "all"), "patchwork")
  expect_s3_class(p, "gg")
})

test_that("internal plot helpers return placeholders for empty data", {
  empty_qc <- structure(
    list(read_counts = list(), quality_scores = list(),
         gc_content = list(), read_lengths = list(), mapping_rate = list()),
    class = "bb_qc"
  )
  expect_s3_class(bambamR:::.plot_quality(empty_qc), "gg")
  expect_s3_class(bambamR:::.plot_gc(empty_qc), "gg")
  expect_s3_class(bambamR:::.plot_length(empty_qc), "gg")
})

test_that(".empty_plot builds a labelled placeholder ggplot", {
  p <- bambamR:::.empty_plot("nothing here")
  expect_s3_class(p, "gg")
})
