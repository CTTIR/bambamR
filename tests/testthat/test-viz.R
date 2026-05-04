test_that("bb_volcano creates ggplot", {
  de <- data.frame(
    gene = paste0("gene", 1:100),
    log2fc = rnorm(100, 0, 2),
    pvalue = 10^(-runif(100, 0, 5)),
    padj = 10^(-runif(100, 0, 4))
  )
  p <- bb_volcano(de)
  expect_s3_class(p, "gg")
})

test_that("bb_volcano validates input", {
  expect_error(bb_volcano(data.frame(a = 1)), "missing required columns")
})

test_that("bb_volcano with label_genes", {
  de <- data.frame(
    gene = paste0("gene", 1:50),
    log2fc = rnorm(50, 0, 2),
    pvalue = 10^(-runif(50, 0, 5)),
    padj = 10^(-runif(50, 0, 4))
  )
  p <- bb_volcano(de, label_genes = c("gene1", "gene5"))
  expect_s3_class(p, "gg")
})

test_that("bb_pca creates ggplot", {
  set.seed(42)
  counts <- matrix(rpois(500, 100), nrow = 50, ncol = 10,
                   dimnames = list(paste0("g", 1:50), paste0("S", 1:10)))
  meta <- data.frame(
    condition = rep(c("A", "B"), each = 5),
    row.names = paste0("S", 1:10)
  )
  p <- bb_pca(counts, meta, color_by = "condition")
  expect_s3_class(p, "gg")
})

test_that("bb_pca validates color_by column", {
  counts <- matrix(rpois(100, 100), nrow = 10,
                   dimnames = list(paste0("g", 1:10), paste0("s", 1:10)))
  meta <- data.frame(grp = rep("A", 10), row.names = paste0("s", 1:10))
  expect_error(bb_pca(counts, meta, color_by = "nonexistent"), "not found")
})

test_that("bb_heatmap creates ggplot", {
  set.seed(42)
  counts <- matrix(rpois(200, 100), nrow = 20, ncol = 10,
                   dimnames = list(paste0("g", 1:20), paste0("s", 1:10)))
  p <- bb_heatmap(counts, n_genes = 10)
  expect_s3_class(p, "gg")
})

test_that("bb_heatmap with DE result", {
  set.seed(42)
  counts <- matrix(rpois(200, 100), nrow = 20, ncol = 10,
                   dimnames = list(paste0("g", 1:20), paste0("s", 1:10)))
  de <- data.frame(
    gene = paste0("g", 1:20),
    log2fc = rnorm(20),
    pvalue = runif(20, 0, 0.1),
    padj = runif(20, 0, 0.1)
  )
  p <- bb_heatmap(counts, de_result = de, n_genes = 10)
  expect_s3_class(p, "gg")
})

test_that("bb_ma_plot creates ggplot", {
  de <- data.frame(
    gene = paste0("gene", 1:100),
    log2fc = rnorm(100, 0, 2),
    pvalue = runif(100, 0, 1),
    padj = runif(100, 0, 1),
    basemean = 10^runif(100, 1, 4)
  )
  p <- bb_ma_plot(de)
  expect_s3_class(p, "gg")
})

test_that("bb_ma_plot requires basemean", {
  de <- data.frame(gene = "A", log2fc = 1, pvalue = 0.01, padj = 0.05)
  expect_error(bb_ma_plot(de), "basemean")
})
