#' PCA Plot
#'
#' Creates a PCA plot from a normalized count matrix with sample metadata.
#'
#' @param counts Numeric matrix. Normalized count matrix (genes x samples).
#' @param metadata A data.frame with sample information. Rownames must match
#'   column names of `counts`.
#' @param color_by Character. Column name in `metadata` to use for coloring
#'   points.
#' @param shape_by Character or NULL. Column name in `metadata` for point
#'   shapes.
#' @param n_genes Integer. Number of top variable genes to use for PCA.
#'   Default `500`.
#' @param label Logical. Whether to label sample points. Default `FALSE`.
#' @param point_size Numeric. Size of points. Default `3`.
#'
#' @return A [ggplot2::ggplot] object.
#'
#' @examples
#' set.seed(42)
#' counts <- matrix(rpois(600, 100), nrow = 100, ncol = 6,
#'   dimnames = list(paste0("gene", 1:100), paste0("S", 1:6)))
#' meta <- data.frame(
#'   condition = rep(c("Control", "Treatment"), each = 3),
#'   row.names = paste0("S", 1:6)
#' )
#' bb_pca(counts, meta, color_by = "condition")
#'
#' @importFrom stats prcomp var
#' @export
bb_pca <- function(counts, metadata, color_by, shape_by = NULL,
                   n_genes = 500L, label = FALSE, point_size = 3) {

  counts <- validate_counts(counts)

  if (!color_by %in% colnames(metadata)) {
    stop("'", color_by, "' not found in metadata columns.", call. = FALSE)
  }

  # Match sample order
  shared <- intersect(colnames(counts), rownames(metadata))
  if (length(shared) == 0L) {
    stop("No matching sample names between counts colnames and metadata rownames.",
         call. = FALSE)
  }
  counts <- counts[, shared, drop = FALSE]
  metadata <- metadata[shared, , drop = FALSE]

  # Select top variable genes
  gene_vars <- apply(counts, 1, stats::var)
  n_genes <- min(n_genes, nrow(counts))
  top_idx <- order(gene_vars, decreasing = TRUE)[seq_len(n_genes)]
  mat <- counts[top_idx, , drop = FALSE]

  # Remove zero-variance rows
  keep <- apply(mat, 1, stats::var) > 0
  mat <- mat[keep, , drop = FALSE]

  if (nrow(mat) < 2L) {
    stop("Not enough variable genes for PCA.", call. = FALSE)
  }

  # PCA
  pca <- stats::prcomp(t(mat), center = TRUE, scale. = TRUE)
  var_explained <- summary(pca)$importance[2, 1:2] * 100

  pca_df <- data.frame(
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    sample_name = rownames(pca$x),
    stringsAsFactors = FALSE
  )
  pca_df[[color_by]] <- metadata[[color_by]]
  if (!is.null(shape_by)) {
    if (!shape_by %in% colnames(metadata)) {
      stop("'", shape_by, "' not found in metadata columns.", call. = FALSE)
    }
    pca_df[[shape_by]] <- metadata[[shape_by]]
  }

  p <- ggplot2::ggplot(pca_df, ggplot2::aes(x = PC1, y = PC2))

  if (!is.null(shape_by)) {
    p <- p + ggplot2::geom_point(
      ggplot2::aes(color = .data[[color_by]], shape = .data[[shape_by]]),
      size = point_size
    )
  } else {
    p <- p + ggplot2::geom_point(
      ggplot2::aes(color = .data[[color_by]]),
      size = point_size
    )
  }

  p <- p +
    ggplot2::labs(
      x = sprintf("PC1 (%.1f%%)", var_explained[1]),
      y = sprintf("PC2 (%.1f%%)", var_explained[2]),
      title = "PCA Plot"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    )

  if (label) {
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_text_repel(
        ggplot2::aes(label = sample_name),
        size = 3
      )
    } else {
      p <- p + ggplot2::geom_text(
        ggplot2::aes(label = sample_name),
        size = 3, vjust = -0.5
      )
    }
  }

  p
}
