#' Export Results to RDS
#'
#' Saves a bambamR result object or any R object as an RDS file.
#'
#' @param result An R object to save (typically a `bb_result` or DE result).
#' @param path Character. File path for the output `.rds` file.
#'
#' @return The `path` invisibly.
#'
#' @examples
#' tmp <- tempfile(fileext = ".rds")
#' counts <- matrix(rpois(100, 50), nrow = 10,
#'   dimnames = list(paste0("g", 1:10), paste0("s", 1:10)))
#' bb_export_rds(counts, tmp)
#' readRDS(tmp)
#'
#' @export
bb_export_rds <- function(result, path) {
  if (missing(path) || !is.character(path) || length(path) != 1L) {
    stop("'path' must be a single character string.", call. = FALSE)
  }
  dir <- dirname(path)
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
  saveRDS(result, file = path)
  message("Saved to: ", path)

  invisible(path)
}


#' Export DE Results to CSV
#'
#' Writes a differential expression result data.frame to a CSV file.
#'
#' @param de_result A data.frame with DE results (columns: gene, log2fc,
#'   pvalue, padj, and optionally basemean).
#' @param path Character. File path for the output `.csv` file.
#'
#' @return The `path` invisibly.
#'
#' @examples
#' de <- data.frame(
#'   gene = paste0("gene", 1:5),
#'   log2fc = rnorm(5),
#'   pvalue = runif(5, 0, 0.1),
#'   padj = runif(5, 0, 0.1)
#' )
#' tmp <- tempfile(fileext = ".csv")
#' bb_export_csv(de, tmp)
#'
#' @export
bb_export_csv <- function(de_result, path) {
  if (missing(path) || !is.character(path) || length(path) != 1L) {
    stop("'path' must be a single character string.", call. = FALSE)
  }
  dir <- dirname(path)
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
  data.table::fwrite(de_result, file = path, sep = ",")
  message("Saved to: ", path)
  invisible(path)
}


#' Export DE Results to TSV
#'
#' Writes a differential expression result data.frame to a TSV file.
#'
#' @param de_result A data.frame with DE results.
#' @param path Character. File path for the output `.tsv` file.
#'
#' @return The `path` invisibly.
#'
#' @examples
#' de <- data.frame(
#'   gene = paste0("gene", 1:5),
#'   log2fc = rnorm(5),
#'   pvalue = runif(5, 0, 0.1),
#'   padj = runif(5, 0, 0.1)
#' )
#' tmp <- tempfile(fileext = ".tsv")
#' bb_export_tsv(de, tmp)
#'
#' @export
bb_export_tsv <- function(de_result, path) {
  if (missing(path) || !is.character(path) || length(path) != 1L) {
    stop("'path' must be a single character string.", call. = FALSE)
  }
  dir <- dirname(path)
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
  data.table::fwrite(de_result, file = path, sep = "\t")
  message("Saved to: ", path)
  invisible(path)
}
