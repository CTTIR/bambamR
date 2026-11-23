# The vectorised QC helpers must reproduce the per-read reference results.

test_that(".compute_gc is vectorised, case-insensitive and unnamed", {
  seqs <- c(a = "GCGC", b = "atgc", c = "ATAT", d = "GgCcN")
  expect_identical(bambamR:::.compute_gc(seqs), c(1, 0.5, 0, 0.8))
  expect_identical(bambamR:::.compute_gc(character(0)), numeric(0))
})

test_that(".compute_qual_summary matches the per-position reference", {
  reference <- function(quality_strings) {
    qual_nums <- lapply(quality_strings, function(q) as.integer(charToRaw(q)) - 33L)
    max_len <- max(lengths(qual_nums))
    vapply(seq_len(max_len), function(pos) {
      mean(vapply(qual_nums, function(q) if (pos <= length(q)) q[pos] else NA_integer_,
                  integer(1)), na.rm = TRUE)
    }, numeric(1))
  }
  set.seed(10)
  symbols <- c("#", "5", "?", "F", "I")
  ragged <- vapply(sample(1:12, 200, replace = TRUE), function(n) {
    paste(sample(symbols, n, replace = TRUE), collapse = "")
  }, character(1))
  res <- bambamR:::.compute_qual_summary(ragged)
  expect_equal(res$position, seq_len(12))
  expect_equal(res$mean_quality, reference(ragged))

  # several accumulation blocks, including reads longer than one block and
  # empty quality strings
  expect_equal(bambamR:::.compute_qual_summary(ragged, block_bases = 5)$mean_quality,
               reference(ragged))
  mixed <- c("", "II#", "", "5", "F?5I#?F", "")
  expect_equal(bambamR:::.compute_qual_summary(mixed, block_bases = 3)$mean_quality,
               reference(mixed))
  many <- rep(c("II#", "5", "F?"), length.out = 50003)
  expect_equal(bambamR:::.compute_qual_summary(many, block_bases = 1e4)$mean_quality,
               reference(many))
  expect_equal(nrow(bambamR:::.compute_qual_summary(c("", ""))), 0L)
})
