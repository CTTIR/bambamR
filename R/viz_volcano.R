#' Volcano Plot
#'
#' Creates a publication-ready volcano plot from differential expression results.
#'
#' @param de_result A data.frame with columns `gene`, `log2fc`, `pvalue`, `padj`.
#'   Typically the output of [bb_deseq2()], [bb_edger()], or [bb_limma_voom()].
#' @param fc_cutoff Numeric. Absolute log2 fold-change cutoff for significance.
#'   Default `1`.
#' @param p_cutoff Numeric. Adjusted p-value cutoff for significance.
#'   Default `0.05`.
#' @param label_genes Character vector. Specific gene names to label on the
#'   plot. If `NULL` (default), the top `n_label` significant genes are labeled.
#' @param n_label Integer. Number of top significant genes to auto-label
#'   when `label_genes` is `NULL`. Default `10`.
#' @param point_size Numeric. Size of points. Default `1`.
#' @param colors Named character vector of length 3 for up, down, and
#'   non-significant colors.
#'
#' @return A [ggplot2::ggplot] object.
#'
#' @examples
#' de <- data.frame(
#'   gene = paste0("gene", 1:200),
#'   log2fc = rnorm(200, 0, 2),
#'   pvalue = 10^(-runif(200, 0, 5)),
#'   padj = 10^(-runif(200, 0, 4))
#' )
#' bb_volcano(de)
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_hline geom_vline labs
#'   scale_color_manual theme_minimal theme element_text
#' @export
bb_volcano <- function(de_result,
                       fc_cutoff = 1,
                       p_cutoff = 0.05,
                       label_genes = NULL,
                       n_label = 10L,
                       point_size = 1,
                       colors = c(up = "#D73027", down = "#4575B4",
                                  ns = "grey70")) {

  validate_de_result(de_result)

  # Remove rows with NA p-values

  df <- de_result[!is.na(de_result$padj), , drop = FALSE]

  # Classify significance
  df$sig <- ifelse(
    df$padj < p_cutoff & df$log2fc > fc_cutoff, "Up",
    ifelse(df$padj < p_cutoff & df$log2fc < -fc_cutoff, "Down", "NS")
  )

  # Label column
  if (!is.null(label_genes)) {
    df$label <- ifelse(df$gene %in% label_genes, df$gene, NA_character_)
  } else {
    sig_genes <- df[df$sig != "NS", , drop = FALSE]
    sig_genes <- sig_genes[order(sig_genes$padj), , drop = FALSE]
    top_genes <- utils::head(sig_genes$gene, n_label)
    df$label <- ifelse(df$gene %in% top_genes, df$gene, NA_character_)
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = log2fc, y = -log10(padj))) +
    ggplot2::geom_point(
      ggplot2::aes(color = sig),
      size = point_size, alpha = 0.7
    ) +
    ggplot2::geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed",
                        color = "grey40") +
    ggplot2::geom_vline(xintercept = c(-fc_cutoff, fc_cutoff),
                        linetype = "dashed", color = "grey40") +
    ggplot2::scale_color_manual(
      values = c("Up" = colors[["up"]], "Down" = colors[["down"]],
                 "NS" = colors[["ns"]]),
      name = "Significance"
    ) +
    ggplot2::labs(
      x = expression(log[2] ~ "Fold Change"),
      y = expression(-log[10] ~ "Adjusted P-value"),
      title = "Volcano Plot"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      legend.position = "right",
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    )

  # Add labels if ggrepel is available
  has_labels <- any(!is.na(df$label))
  if (has_labels && requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p + ggrepel::geom_text_repel(
      ggplot2::aes(label = label),
      size = 3, max.overlaps = 20,
      segment.color = "grey50"
    )
  } else if (has_labels) {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = label),
      size = 3, vjust = -0.5, hjust = 0.5, check_overlap = TRUE
    )
  }

  p
}
