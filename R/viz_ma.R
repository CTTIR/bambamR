#' MA Plot
#'
#' Creates an MA plot (log2 fold-change vs. mean expression) from
#' differential expression results.
#'
#' @param de_result A data.frame with columns `gene`, `log2fc`, `pvalue`,
#'   `padj`, and `basemean`.
#' @param p_cutoff Numeric. Adjusted p-value cutoff for coloring significant
#'   genes. Default `0.05`.
#' @param point_size Numeric. Size of points. Default `1`.
#' @param colors Named character vector with `"sig"` and `"ns"` colors.
#'
#' @return A [ggplot2::ggplot] object.
#'
#' @examples
#' de <- data.frame(
#'   gene = paste0("gene", 1:200),
#'   log2fc = rnorm(200, 0, 2),
#'   pvalue = 10^(-runif(200, 0, 5)),
#'   padj = 10^(-runif(200, 0, 4)),
#'   basemean = 10^runif(200, 1, 4)
#' )
#' bb_ma_plot(de)
#'
#' @export
bb_ma_plot <- function(de_result, p_cutoff = 0.05, point_size = 1,
                       colors = c(sig = "#D73027", ns = "grey70")) {

  validate_de_result(de_result)

  if (!"basemean" %in% colnames(de_result)) {
    stop("'de_result' must contain a 'basemean' column for MA plots.",
         call. = FALSE)
  }

  df <- de_result[!is.na(de_result$padj) & !is.na(de_result$basemean), ,
                  drop = FALSE]

  df$sig <- ifelse(df$padj < p_cutoff, "Significant", "NS")

  p <- ggplot2::ggplot(df, ggplot2::aes(x = log10(basemean), y = log2fc)) +
    ggplot2::geom_point(
      ggplot2::aes(color = sig),
      size = point_size, alpha = 0.6
    ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey30") +
    ggplot2::scale_color_manual(
      values = c("Significant" = colors[["sig"]], "NS" = colors[["ns"]]),
      name = NULL
    ) +
    ggplot2::labs(
      x = expression(log[10] ~ "Mean Expression"),
      y = expression(log[2] ~ "Fold Change"),
      title = "MA Plot"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      legend.position = "right"
    )

  p
}
