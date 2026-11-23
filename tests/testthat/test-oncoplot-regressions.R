# Regression tests for oncoplot layout, palettes and argument validation.

onco_data <- function() {
  set.seed(11)
  data.frame(
    sample = sample(paste0("S", 1:12), 60, replace = TRUE),
    gene = sample(c("TP53", "KRAS", "PTEN", "APC", "BRAF"), 60, replace = TRUE),
    mutation_type = sample(c("Missense_Mutation", "Nonsense_Mutation",
                             "Splice_Site"), 60, replace = TRUE),
    stringsAsFactors = FALSE
  )
}

test_that("bb_oncoplot accepts a custom palette without an 'Other' entry", {
  dat <- data.frame(
    sample = c("S1", "S2", "S3"),
    gene = c("TP53", "KRAS", "PTEN"),
    mutation_type = c("Unlisted_Type", "Missense_Mutation", "Splice_Site"),
    stringsAsFactors = FALSE
  )
  p <- bb_oncoplot(dat, show_barplot = FALSE,
                   mutation_colors = c(Missense_Mutation = "#0072B2",
                                       Nonsense_Mutation = "#D55E00",
                                       Splice_Site = "#E69F00"))
  fill <- ggplot2::ggplot_build(p)$plot$scales$get_scales("fill")
  expect_equal(unname(fill$map(c("Unlisted_Type", "Splice_Site"))),
               c("#999999", "#E69F00"))
  expect_error(bb_oncoplot(dat, mutation_colors = c("red", "blue")),
               "named character vector")
})

test_that("bb_oncoplot shows duplicated 'genes' once and validates gene selection", {
  dat <- onco_data()
  p <- bb_oncoplot(dat, genes = c("TP53", "KRAS", "TP53"), show_barplot = FALSE)
  expect_equal(sort(levels(p$data$gene)), c("KRAS", "TP53"))
  # a factor of gene names keeps working, as before the validation was added
  p_factor <- bb_oncoplot(dat, genes = factor(c("TP53", "KRAS")), show_barplot = FALSE)
  expect_equal(sort(levels(p_factor$data$gene)), c("KRAS", "TP53"))
  expect_error(bb_oncoplot(dat, genes = c("TP53", NA)), "'genes' must be")
  expect_error(bb_oncoplot(dat, n_genes = 0), "'n_genes' must be a positive")
  expect_error(bb_oncoplot(dat, n_genes = -2), "'n_genes' must be a positive")
  expect_error(bb_oncoplot(dat, n_genes = NA_integer_), "'n_genes' must be")
  all_genes <- bb_oncoplot(dat, n_genes = Inf, show_barplot = FALSE)
  expect_equal(nlevels(all_genes$data$gene), 5L)
})

test_that("bb_oncoplot warns when no annotation row matches a sample", {
  dat <- onco_data()
  anno <- data.frame(stage = c("I", "II"), row.names = c("other1", "other2"))
  expect_warning(
    p <- bb_oncoplot(dat, annotation_df = anno, show_barplot = FALSE),
    "annotation tracks are omitted"
  )
  expect_s3_class(p, "gg")
  expect_error(bb_oncoplot(dat, annotation_df = "stage"), "'annotation_df' must be")
})

test_that("bb_oncoplot collects legends below the whole composition", {
  skip_if_not_installed("patchwork")
  p <- bb_oncoplot(onco_data())
  expect_s3_class(p, "patchwork")
  expect_identical(p$patches$layout$guides, "collect")
  # the top barplot title is broken over two lines so it fits the short panel
  y_titles <- vapply(seq_along(p$patches$plots) + 1L, function(i) {
    title <- ggplot2::ggplot_build(p[[i]])$plot$labels$y
    if (is.null(title)) NA_character_ else as.character(title)
  }, character(1))
  expect_true("Mutated\ngenes" %in% y_titles)
})

test_that("the top barplot of gene counts uses whole-number axis breaks", {
  skip_if_not_installed("patchwork")
  p <- bb_oncoplot(bb_example_mutations()$mutations, n_genes = 20)
  # patchwork orders the panels by area name, so find the top barplot by title
  built <- lapply(seq_len(length(p$patches$plots) + 1L),
                  function(i) ggplot2::ggplot_build(p[[i]]))
  is_top <- vapply(built, function(b) identical(b$plot$labels$y, "Mutated\ngenes"),
                   logical(1))
  expect_equal(sum(is_top), 1L)
  breaks <- built[[which(is_top)]]$layout$panel_params[[1]]$y$breaks
  breaks <- breaks[!is.na(breaks)]
  expect_true(is.numeric(breaks))
  expect_gte(length(breaks), 2L)
  expect_equal(breaks, round(breaks))
})

test_that("cluster ordering of oncoplot samples matches the row-wise reference", {
  dat <- bambamR:::.handle_multi_hit(onco_data())
  genes <- unique(dat$gene)
  samples <- unique(dat$sample)
  mat <- matrix(0L, nrow = length(genes), ncol = length(samples),
                dimnames = list(genes, samples))
  for (i in seq_len(nrow(dat))) mat[dat$gene[i], dat$sample[i]] <- 1L
  reference <- colnames(mat)[stats::hclust(stats::dist(t(mat), method = "binary"))$order]
  expect_identical(bambamR:::.sort_onco_samples(dat, genes, "cluster"), reference)
})

test_that("annotation tracks without barplots also collect their legends", {
  skip_if_not_installed("patchwork")
  dat <- onco_data()
  anno <- data.frame(stage = rep(c("I", "II"), 6), row.names = paste0("S", 1:12))
  p <- bb_oncoplot(dat, annotation_df = anno, show_barplot = FALSE)
  expect_s3_class(p, "patchwork")
  expect_identical(p$patches$layout$guides, "collect")
})
