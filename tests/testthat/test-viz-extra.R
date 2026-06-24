# Extra visualization tests covering uncovered branches of R/viz_pca.R,
# R/viz_heatmap.R and R/viz_volcano.R.

test_that("bb_pca supports shape_by and labelling", {
  set.seed(21)
  counts <- matrix(rpois(600, 100), nrow = 100, ncol = 6,
                   dimnames = list(paste0("g", 1:100), paste0("S", 1:6)))
  meta <- data.frame(
    condition = rep(c("A", "B"), each = 3),
    batch = rep(c("x", "y", "z"), 2),
    row.names = paste0("S", 1:6)
  )
  p <- bb_pca(counts, meta, color_by = "condition", shape_by = "batch",
              label = TRUE)
  expect_s3_class(p, "gg")
})

test_that("bb_pca errors on an unknown shape_by column", {
  counts <- matrix(rpois(200, 100), nrow = 20, ncol = 10,
                   dimnames = list(paste0("g", 1:20), paste0("s", 1:10)))
  meta <- data.frame(condition = rep(c("A", "B"), 5),
                     row.names = paste0("s", 1:10))
  expect_error(
    bb_pca(counts, meta, color_by = "condition", shape_by = "nope"),
    "not found"
  )
})

test_that("bb_pca errors when no sample names overlap", {
  counts <- matrix(rpois(200, 100), nrow = 20, ncol = 10,
                   dimnames = list(paste0("g", 1:20), paste0("s", 1:10)))
  meta <- data.frame(condition = rep("A", 10),
                     row.names = paste0("x", 1:10))
  expect_error(bb_pca(counts, meta, color_by = "condition"), "No matching")
})

test_that("bb_heatmap honours scale = 'none' and 'column'", {
  set.seed(22)
  counts <- matrix(rpois(200, 100), nrow = 20, ncol = 10,
                   dimnames = list(paste0("g", 1:20), paste0("s", 1:10)))
  expect_s3_class(bb_heatmap(counts, n_genes = 10, scale = "none"), "gg")
  expect_s3_class(bb_heatmap(counts, n_genes = 10, scale = "column"), "gg")
})

test_that("bb_heatmap can disable clustering", {
  set.seed(23)
  counts <- matrix(rpois(200, 100), nrow = 20, ncol = 10,
                   dimnames = list(paste0("g", 1:20), paste0("s", 1:10)))
  p <- bb_heatmap(counts, n_genes = 10, cluster_rows = FALSE,
                  cluster_cols = FALSE)
  expect_s3_class(p, "gg")
})

test_that("bb_heatmap errors when DE genes do not match counts", {
  set.seed(24)
  counts <- matrix(rpois(200, 100), nrow = 20, ncol = 10,
                   dimnames = list(paste0("g", 1:20), paste0("s", 1:10)))
  de <- data.frame(gene = paste0("other", 1:5),
                   log2fc = rnorm(5), pvalue = runif(5), padj = runif(5))
  expect_error(bb_heatmap(counts, de_result = de), "No matching genes")
})

test_that("bb_volcano labels fall back to geom_text without ggrepel", {
  de <- data.frame(
    gene = paste0("gene", 1:50),
    log2fc = c(rep(5, 5), rnorm(45)),
    pvalue = c(rep(1e-6, 5), runif(45, 0.1, 1)),
    padj = c(rep(1e-5, 5), runif(45, 0.1, 1))
  )
  local_mocked_bindings(
    requireNamespace = function(pkg, ...) {
      if (identical(pkg, "ggrepel")) FALSE else TRUE
    },
    .package = "base"
  )
  p <- bb_volcano(de, label_genes = c("gene1", "gene2"))
  expect_s3_class(p, "gg")
})

test_that("bb_volcano auto-labels top significant genes", {
  de <- data.frame(
    gene = paste0("gene", 1:50),
    log2fc = c(rep(5, 5), rep(-5, 5), rnorm(40)),
    pvalue = c(rep(1e-8, 10), runif(40, 0.5, 1)),
    padj = c(rep(1e-7, 10), runif(40, 0.5, 1))
  )
  p <- bb_volcano(de)
  expect_s3_class(p, "gg")
})
