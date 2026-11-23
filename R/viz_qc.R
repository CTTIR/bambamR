#' QC Visualization Panel
#'
#' Creates a panel of QC plots from a `bb_qc` object, including read quality
#' per position, GC content distribution, and read length distribution.
#'
#' @param qc_object A `bb_qc` object from [bb_qc()].
#' @param which Character. Which plot to generate: `"quality"`, `"gc"`,
#'   `"length"`, or `"all"`. Default `"all"`.
#'
#' @return A [ggplot2::ggplot] object (or a patchwork composition if
#'   `which = "all"` and patchwork is available). In the composition the file
#'   legend is drawn once, below the last panel; files keep the same color in
#'   every panel.
#'
#' @examples
#' \donttest{
#' tmp <- tempfile(fileext = ".fastq")
#' writeLines(c(
#'   "@r1", "ACGTACGT", "+", "IIIIIIII",
#'   "@r2", "GCGCGCGC", "+", "HHHHHHHH",
#'   "@r3", "ATATGCGC", "+", "GGGGFFFF"
#' ), tmp)
#' qc <- bb_qc(fastq_path = tmp)
#' bb_plot_qc(qc)
#' }
#'
#' @export
bb_plot_qc <- function(qc_object, which = c("all", "quality", "gc", "length")) {

  if (!methods::is(qc_object, "bb_qc")) {
    stop("'qc_object' must be a bb_qc object.", call. = FALSE)
  }

  which <- match.arg(which)

  plots <- list()

  if (which %in% c("all", "quality")) {
    plots$quality <- .plot_quality(qc_object)
  }
  if (which %in% c("all", "gc")) {
    plots$gc <- .plot_gc(qc_object)
  }
  if (which %in% c("all", "length")) {
    plots$length <- .plot_length(qc_object)
  }

  if (length(plots) == 1L) {
    return(plots[[1]])
  }

  if (requireNamespace("patchwork", quietly = TRUE)) {
    # Files have the same color in every panel, so one legend (below the last
    # panel) is enough; three legends squash the panels at common figure sizes.
    upper <- seq_len(length(plots) - 1L)
    plots[upper] <- lapply(plots[upper], function(p) {
      p + ggplot2::theme(legend.position = "none")
    })
    combined <- Reduce(`+`, plots) +
      patchwork::plot_layout(ncol = 1)
    return(combined)
  }

  # Fallback: return quality plot only with message
  message("Install 'patchwork' for combined QC panels. Returning quality plot.")
  plots[[1]]
}


#' Per-position quality plot
#' @noRd
.plot_quality <- function(qc_object) {
  qual_list <- qc_object$quality_scores
  if (length(qual_list) == 0L) {
    return(.empty_plot("No quality data available"))
  }

  # Combine all files
  dfs <- lapply(names(qual_list), function(nm) {
    df <- qual_list[[nm]]
    if (is.null(df)) return(NULL)
    df$file <- nm
    df
  })
  dfs <- dfs[!vapply(dfs, is.null, logical(1))]
  if (length(dfs) == 0L) return(.empty_plot("No quality data available"))

  df_all <- do.call(rbind, dfs)

  ggplot2::ggplot(df_all, ggplot2::aes(x = position, y = mean_quality)) +
    ggplot2::geom_line(ggplot2::aes(color = file), linewidth = 0.8) +
    # Thresholds are labelled and differ in line type as well as in color, so
    # they do not rely on a red-green contrast. Colors: vermillion (#D55E00)
    # and bluish green (#009E73) of the color-blind-safe palette in Wong
    # (2011), Nat Methods 8(6):441, doi:10.1038/nmeth.1618, Figure 2, p. 441.
    ggplot2::geom_hline(yintercept = 20, linetype = "dashed", color = "#D55E00") +
    ggplot2::geom_hline(yintercept = 30, linetype = "dotdash", color = "#009E73") +
    ggplot2::annotate("text", x = Inf, y = c(20, 30), label = c("Q20", "Q30"),
                      hjust = 1.1, vjust = -0.4, size = 3,
                      color = c("#D55E00", "#009E73")) +
    ggplot2::labs(
      x = "Position in Read",
      y = "Mean Phred\nQuality",
      title = "Per-Position Quality",
      color = "File"
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      legend.position = "bottom"
    )
}


#' GC content distribution plot
#' @noRd
.plot_gc <- function(qc_object) {
  gc_list <- qc_object$gc_content
  if (length(gc_list) == 0L) {
    return(.empty_plot("No GC content data available"))
  }

  dfs <- lapply(names(gc_list), function(nm) {
    gc <- gc_list[[nm]]
    if (is.null(gc) || length(gc) == 0L) return(NULL)
    data.frame(gc_content = gc, file = nm, stringsAsFactors = FALSE)
  })
  dfs <- dfs[!vapply(dfs, is.null, logical(1))]
  if (length(dfs) == 0L) return(.empty_plot("No GC content data available"))

  df_all <- do.call(rbind, dfs)

  ggplot2::ggplot(df_all, ggplot2::aes(x = gc_content)) +
    ggplot2::geom_density(ggplot2::aes(fill = file), alpha = 0.4) +
    ggplot2::labs(
      x = "GC Content",
      y = "Density",
      title = "GC Content Distribution",
      fill = "File"
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      legend.position = "bottom"
    )
}


#' Read length distribution plot
#' @noRd
.plot_length <- function(qc_object) {
  len_list <- qc_object$read_lengths
  if (length(len_list) == 0L) {
    return(.empty_plot("No read length data available"))
  }

  dfs <- lapply(names(len_list), function(nm) {
    lens <- len_list[[nm]]
    if (is.null(lens) || length(lens) == 0L) return(NULL)
    data.frame(read_length = lens, file = nm, stringsAsFactors = FALSE)
  })
  dfs <- dfs[!vapply(dfs, is.null, logical(1))]
  if (length(dfs) == 0L) return(.empty_plot("No read length data available"))

  df_all <- do.call(rbind, dfs)

  # Read lengths are whole numbers: one bin per length while the range is
  # narrow (also when every read has the same length), 50 bins otherwise.
  histogram <- if (diff(range(df_all$read_length)) < 50) {
    ggplot2::geom_histogram(ggplot2::aes(fill = file), binwidth = 1,
                            alpha = 0.6, position = "identity")
  } else {
    ggplot2::geom_histogram(ggplot2::aes(fill = file), bins = 50,
                            alpha = 0.6, position = "identity")
  }
  whole_breaks <- function(limits) {
    breaks <- unique(round(pretty(limits)))
    breaks[breaks >= limits[1] & breaks <= limits[2]]
  }

  ggplot2::ggplot(df_all, ggplot2::aes(x = read_length)) +
    histogram +
    ggplot2::scale_x_continuous(
      breaks = whole_breaks,
      expand = ggplot2::expansion(mult = 0.05, add = 0.5)
    ) +
    ggplot2::labs(
      x = "Read Length",
      y = "Count",
      title = "Read Length Distribution",
      fill = "File"
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      legend.position = "bottom"
    )
}


#' Empty placeholder plot
#' @noRd
.empty_plot <- function(msg) {
  ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0.5, y = 0.5, label = msg, size = 5) +
    ggplot2::theme_void()
}
