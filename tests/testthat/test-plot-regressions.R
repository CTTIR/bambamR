# Regression tests for heatmap, volcano, MA, PCA and QC plot defects.

de_data <- function() {
  set.seed(3)
  data.frame(
    gene = paste0("g", 1:60),
    log2fc = c(rep(3, 5), rep(-3, 5), rnorm(50, 0, 0.3)),
    pvalue = c(rep(1e-6, 10), runif(50, 0.2, 1)),
    padj = c(seq(1e-5, 1e-4, length.out = 10), runif(50, 0.3, 1)),
    basemean = 10^runif(60, 1, 4),
    stringsAsFactors = FALSE
  )
}

# Heatmap -------------------------------------------------------------

test_that("bb_heatmap centres the Z-score color scale at zero", {
  set.seed(5)
  counts <- matrix(rpois(120, 80), nrow = 20,
                   dimnames = list(paste0("g", 1:20), paste0("S", 1:6)))
  counts[1, 6] <- 5000
  p <- bb_heatmap(counts, n_genes = 20)
  limits <- p$scales$get_scales("fill")$limits
  expect_equal(limits[1], -limits[2])
  expect_equal(limits[2], max(abs(p$data$value)))
  expect_null(bb_heatmap(counts, scale = "none")$scales$get_scales("fill")$limits)
})

test_that("bb_heatmap shows a single sample as zero instead of missing values", {
  counts <- matrix(c(5, 10, 20), ncol = 1,
                   dimnames = list(c("g1", "g2", "g3"), "S1"))
  p <- bb_heatmap(counts)
  expect_false(anyNA(p$data$value))
  expect_true(all(p$data$value == 0))
})

test_that("bb_heatmap validates n_genes and labels samples once with annotations", {
  set.seed(6)
  counts <- matrix(rpois(60, 80), nrow = 10,
                   dimnames = list(paste0("g", 1:10), paste0("S", 1:6)))
  expect_error(bb_heatmap(counts, n_genes = 0), "'n_genes' must be a positive")
  expect_error(bb_heatmap(counts, n_genes = "5"), "'n_genes' must be a positive")
  expect_error(bb_heatmap(counts, n_genes = 2.5), "'n_genes' must be a positive")
  # Inf keeps selecting every gene
  expect_equal(nlevels(bb_heatmap(counts, n_genes = Inf)$data$row_gene), 10L)
  skip_if_not_installed("patchwork")
  anno <- data.frame(group = rep(c("a", "b"), 3), row.names = colnames(counts))
  p <- bb_heatmap(counts, annotation_col = anno)
  expect_s3_class(p[[1]]$theme$axis.text.x, "element_blank")
  expect_s3_class(p[[2]]$theme$axis.text.x, "element_text")
})

# Volcano and MA ------------------------------------------------------

test_that("bb_volcano draws without a missing-value warning for unlabelled genes", {
  de <- de_data()
  p <- bb_volcano(de, n_label = 4)
  label_layer <- p$layers[[length(p$layers)]]
  expect_equal(nrow(label_layer$data), 4L)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_warning(ggplot2::ggplotGrob(p))
})

test_that("bb_volcano validates colors, n_label and drops unplottable rows", {
  de <- de_data()
  n_layers <- length(bb_volcano(de, n_label = 0)$layers)
  expect_equal(n_layers, 3L)
  expect_s3_class(bb_volcano(de, colors = c("red", "blue", "grey")), "gg")
  expect_error(bb_volcano(de, colors = c(high = "red", low = "blue", ns = "grey")),
               "named \"up\", \"down\", \"ns\"")
  expect_error(bb_volcano(de, n_label = -1), "'n_label' must be")
  inf_labels <- bb_volcano(de, n_label = Inf)
  expect_equal(nrow(inf_labels$layers[[length(inf_labels$layers)]]$data), 10L)
  de$log2fc[1] <- NA
  expect_equal(nrow(bb_volcano(de)$data), nrow(de) - 1L)
})

test_that("bb_ma_plot reads unnamed colors in order and rejects wrong names", {
  de <- de_data()
  p <- bb_ma_plot(de, colors = c("black", "grey80"))
  colour <- ggplot2::ggplot_build(p)$plot$scales$get_scales("colour")
  expect_equal(unname(colour$map(c("Significant", "NS"))), c("black", "grey80"))
  expect_error(bb_ma_plot(de, colors = c(a = "black", b = "grey")),
               "named \"sig\", \"ns\"")
})

# PCA -----------------------------------------------------------------

test_that("bb_pca validates metadata, color_by and n_genes", {
  set.seed(8)
  counts <- matrix(rpois(60, 80), nrow = 10,
                   dimnames = list(paste0("g", 1:10), paste0("S", 1:6)))
  meta <- data.frame(group = rep(c("a", "b"), 3), row.names = colnames(counts))
  expect_error(bb_pca(counts, as.matrix(meta), color_by = "group"),
               "'metadata' must be a data.frame")
  expect_error(bb_pca(counts, meta, color_by = c("group", "group")),
               "'color_by' must be a single column name")
  expect_error(bb_pca(counts, meta, color_by = "group", n_genes = 0),
               "'n_genes' must be a positive")
})

# QC plots ------------------------------------------------------------

test_that("read-length histogram uses whole-number bins and breaks", {
  qc <- structure(list(read_counts = list(a = 3), quality_scores = list(),
                       gc_content = list(), read_lengths = list(a = c(8L, 8L, 8L)),
                       mapping_rate = list()), class = "bb_qc")
  p <- bb_plot_qc(qc, which = "length")
  built <- ggplot2::ggplot_build(p)
  expect_equal(built$data[[1]]$xmax - built$data[[1]]$xmin, 1)
  breaks <- built$layout$panel_params[[1]]$x$breaks
  breaks <- breaks[!is.na(breaks)]
  expect_true(length(breaks) >= 1L)
  expect_equal(breaks, round(breaks))

  qc$read_lengths$a <- c(30L, 90L, 151L)
  wide <- ggplot2::ggplot_build(bb_plot_qc(qc, which = "length"))
  expect_equal(nrow(wide$data[[1]]), 50L)
})

test_that("quality thresholds are labelled in the per-position plot", {
  qc <- structure(list(
    read_counts = list(a = 2),
    quality_scores = list(a = data.frame(position = 1:4, mean_quality = c(38, 37, 36, 35))),
    gc_content = list(), read_lengths = list(), mapping_rate = list()
  ), class = "bb_qc")
  p <- bb_plot_qc(qc, which = "quality")
  is_text <- vapply(p$layers, function(l) inherits(l$geom, "GeomText"), logical(1))
  expect_equal(sum(is_text), 1L)
  expect_setequal(p$layers[[which(is_text)]]$aes_params$label, c("Q20", "Q30"))
})

test_that("the combined QC panel draws the file legend once", {
  skip_if_not_installed("patchwork")
  qc <- structure(list(
    read_counts = list(a = 3),
    quality_scores = list(a = data.frame(position = 1:4, mean_quality = c(38, 37, 36, 35))),
    gc_content = list(a = c(0.4, 0.5, 0.6)),
    read_lengths = list(a = c(4L, 4L, 4L)),
    mapping_rate = list()
  ), class = "bb_qc")
  p <- bb_plot_qc(qc)
  expect_s3_class(p, "patchwork")
  positions <- vapply(1:3, function(i) {
    position <- p[[i]]$theme$legend.position
    if (is.null(position)) NA_character_ else as.character(position)
  }, character(1))
  expect_equal(sum(positions == "none", na.rm = TRUE), 2L)
  expect_equal(sum(positions == "bottom", na.rm = TRUE), 1L)
})
