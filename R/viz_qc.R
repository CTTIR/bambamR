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
#'   `which = "all"` and patchwork is available).
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
    ggplot2::geom_hline(yintercept = 20, linetype = "dashed", color = "red") +
    ggplot2::geom_hline(yintercept = 30, linetype = "dashed", color = "green3") +
    ggplot2::labs(
      x = "Position in Read",
      y = "Mean Phred Quality",
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

  ggplot2::ggplot(df_all, ggplot2::aes(x = read_length)) +
    ggplot2::geom_histogram(ggplot2::aes(fill = file),
                            bins = 50, alpha = 0.6,
                            position = "identity") +
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
