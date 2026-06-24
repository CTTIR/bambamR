# Extra tests for R/oncoplot.R covering sort_by = "cluster", title,
# unknown-type remapping, the no-matching-genes error, and internal helpers.

make_mut <- function(seed = 31, n = 60) {
  set.seed(seed)
  data.frame(
    sample = sample(paste0("S", 1:12), n, replace = TRUE),
    gene = sample(c("TP53", "KRAS", "PIK3CA", "PTEN", "APC", "EGFR"),
                  n, replace = TRUE),
    mutation_type = sample(c("Missense_Mutation", "Nonsense_Mutation",
                             "Splice_Site"), n, replace = TRUE),
    stringsAsFactors = FALSE
  )
}

test_that("bb_oncoplot supports cluster sorting", {
  p <- bb_oncoplot(make_mut(), n_genes = 5, sort_by = "cluster",
                   show_barplot = FALSE)
  expect_true(inherits(p, "gg") || inherits(p, "patchwork"))
})

test_that("bb_oncoplot adds a title in the fallback (no barplot) path", {
  p <- bb_oncoplot(make_mut(), n_genes = 4, show_barplot = FALSE,
                   title = "My Cohort")
  expect_s3_class(p, "gg")
})

test_that("bb_oncoplot remaps unknown mutation types to Other", {
  dat <- data.frame(
    sample = c("S1", "S2", "S3"),
    gene = c("TP53", "KRAS", "PTEN"),
    mutation_type = c("Weird_Type", "Missense_Mutation", "Another_Weird"),
    stringsAsFactors = FALSE
  )
  p <- bb_oncoplot(dat, n_genes = 3, show_barplot = FALSE)
  expect_s3_class(p, "gg")
})

test_that("bb_oncoplot errors when specified genes have no mutations", {
  expect_error(
    bb_oncoplot(make_mut(), genes = c("NONEXISTENT1", "NONEXISTENT2")),
    "No mutations found"
  )
})

test_that("bb_oncoplot rejects an unknown sort_by", {
  expect_error(bb_oncoplot(make_mut(), sort_by = "alpha"), "should be one of")
})

test_that(".standardize_onco_input drops empty sample/gene rows", {
  dat <- data.frame(
    sample = c("S1", "", "S3"),
    gene = c("TP53", "KRAS", ""),
    mutation_type = c("Missense_Mutation", "Missense_Mutation",
                      "Nonsense_Mutation"),
    stringsAsFactors = FALSE
  )
  res <- bambamR:::.standardize_onco_input(dat)
  expect_equal(nrow(res), 1L)
  expect_equal(res$sample, "S1")
})

test_that(".handle_multi_hit collapses repeated sample-gene pairs", {
  df <- data.frame(
    sample = c("S1", "S1", "S2"),
    gene = c("TP53", "TP53", "KRAS"),
    mutation_type = c("Missense_Mutation", "Nonsense_Mutation",
                      "Missense_Mutation"),
    stringsAsFactors = FALSE
  )
  res <- bambamR:::.handle_multi_hit(df)
  multi <- res[res$sample == "S1" & res$gene == "TP53", ]
  expect_equal(multi$mutation_type, "Multi_Hit")
})

test_that(".melt_annotation reshapes annotation tracks to long format", {
  anno <- data.frame(Stage = c("I", "II"), Sex = c("M", "F"),
                     row.names = c("S1", "S2"))
  long <- bambamR:::.melt_annotation(anno, c("S1", "S2"), c("S1", "S2"))
  expect_s3_class(long, "data.frame")
  expect_equal(nrow(long), 4L)
  expect_true(all(c("sample", "variable", "annotation_value") %in%
                    colnames(long)))
})
