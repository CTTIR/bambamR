#' Heatmap of Top Differentially Expressed Genes
#'
#' Creates a heatmap from a normalized count matrix, optionally highlighting
#' the top DE genes.
#'
#' @param counts Numeric matrix. Normalized count matrix (genes x samples).
#' @param de_result A data.frame with DE results. If provided, the top
#'   `n_genes` by adjusted p-value are shown.
#' @param n_genes Integer. Number of top genes to display. Must be a positive
#'   whole number, or `Inf` for all genes. Default `50`.
#' @param annotation_col A data.frame for column (sample) annotations.
#'   Rownames must match column names of `counts`. Annotation tracks require
#'   the optional `patchwork` package.
#' @param scale Character. Scale rows (`"row"`), columns (`"column"`),
#'   or neither (`"none"`). Default `"row"`.
#' @param cluster_rows Logical. Cluster rows. Default `TRUE`.
#' @param cluster_cols Logical. Cluster columns. Default `TRUE`.
#' @param color_palette Character vector. Colours for the heatmap gradient.
#'   The default is a blue-white-red diverging palette.
#'
#' @details
#' With `scale = "row"` or `scale = "column"` the colour scale is symmetric
#' around zero, so the middle colour of `color_palette` always marks a Z-score
#' of zero. Rows or columns without variance are shown as zero.
#'
#' @return A [ggplot2::ggplot] object (a patchwork composition when
#'   `annotation_col` is supplied).
#'
#' @examples
#' set.seed(42)
#' counts <- matrix(rpois(500, 100), nrow = 50, ncol = 10,
#'   dimnames = list(paste0("gene", 1:50), paste0("S", 1:10)))
#' bb_heatmap(counts, n_genes = 20)
#'
#' @importFrom stats hclust dist
#' @importFrom grDevices colorRampPalette
#' @export
bb_heatmap <- function(counts, de_result = NULL, n_genes = 50L,
                       annotation_col = NULL,
                       scale = c("row", "column", "none"),
                       cluster_rows = TRUE, cluster_cols = TRUE,
                       color_palette = NULL) {

  scale <- match.arg(scale)
  counts <- validate_counts(counts)
  .validate_n(n_genes, "n_genes")

  # Select genes
  if (!is.null(de_result)) {
    validate_de_result(de_result)
    de_sorted <- de_result[order(de_result$padj), , drop = FALSE]
    top_genes <- utils::head(de_sorted$gene, n_genes)
    top_genes <- intersect(top_genes, rownames(counts))
    if (length(top_genes) == 0L) {
      stop("No matching genes between DE result and count matrix.", call. = FALSE)
    }
    mat <- counts[top_genes, , drop = FALSE]
  } else {
    # Use top variable genes
    gene_vars <- apply(counts, 1, stats::var)
    n_genes <- min(n_genes, nrow(counts))
    top_idx <- order(gene_vars, decreasing = TRUE)[seq_len(n_genes)]
    mat <- counts[top_idx, , drop = FALSE]
  }

  # Scale
  if (scale == "row") {
    mat <- t(scale(t(mat)))
  } else if (scale == "column") {
    mat <- scale(mat)
  }

  # Zero-variance rows give NaN and a single sample gives NA: show both as 0
  mat[!is.finite(mat)] <- 0

  # Clustering order
  if (cluster_rows && nrow(mat) > 2L) {
    row_order <- stats::hclust(stats::dist(mat))$order
    mat <- mat[row_order, , drop = FALSE]
  }
  if (cluster_cols && ncol(mat) > 2L) {
    col_order <- stats::hclust(stats::dist(t(mat)))$order
    mat <- mat[, col_order, drop = FALSE]
  }

  # Convert to long format for ggplot
  row_genes <- rownames(mat)
  col_samples <- colnames(mat)

  df <- data.frame(
    row_gene = rep(row_genes, times = ncol(mat)),
    col_sample = rep(col_samples, each = nrow(mat)),
    value = as.vector(mat),
    stringsAsFactors = FALSE
  )

  # Preserve ordering

  df$row_gene <- factor(df$row_gene, levels = rev(row_genes))
  df$col_sample <- factor(df$col_sample, levels = col_samples)

  # Color palette
  if (is.null(color_palette)) {
    color_palette <- c("#2166AC", "#F7F7F7", "#B2182B")
  }
  fill_gradient <- grDevices::colorRampPalette(color_palette)(100)
  # Z-scores: symmetric limits keep the palette midpoint at zero.
  fill_limits <- NULL
  if (scale != "none") {
    bound <- max(abs(df$value))
    if (is.finite(bound) && bound > 0) fill_limits <- c(-bound, bound)
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = col_sample, y = row_gene,
                                        fill = value)) +
    ggplot2::geom_tile(color = NA) +
    ggplot2::scale_fill_gradientn(colors = fill_gradient, limits = fill_limits,
                                  name = if (scale == "none") "Expression" else "Z-score") +
    ggplot2::labs(x = NULL, y = NULL, title = "Heatmap") +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = ggplot2::element_text(size = 6),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      panel.grid = ggplot2::element_blank()
    )

  if (!is.null(annotation_col)) {
    check_pkg("patchwork", reason = "for heatmap annotation tracks")
    annotation_col <- .match_metadata(counts, annotation_col)
    anno <- .melt_annotation(annotation_col, col_samples, col_samples)
    p_anno <- ggplot2::ggplot(anno, ggplot2::aes(x = sample, y = variable,
                                                 fill = annotation_value)) +
      ggplot2::geom_tile() +
      ggplot2::scale_x_discrete(limits = col_samples, drop = FALSE) +
      ggplot2::labs(x = NULL, y = NULL, fill = "Annotation") +
      ggplot2::theme_minimal(base_size = 10) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8),
        panel.grid = ggplot2::element_blank()
      )
    # Sample labels are drawn once, below the annotation track.
    p <- p + ggplot2::theme(axis.text.x = ggplot2::element_blank())
    return(p / p_anno + patchwork::plot_layout(heights = c(4, 1)))
  }
  p
}
