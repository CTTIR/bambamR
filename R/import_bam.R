#' Read BAM Files
#'
#' Reads BAM files using Rsamtools (if available) or the system `samtools`
#' command as fallback.
#'
#' @param path Character. Path to a BAM file.
#' @param index Logical. Whether to use/create a BAM index. Default `TRUE`.
#' @param what Character vector. Fields to extract. Default includes basic
#'   alignment fields.
#'
#' @return A data.frame with alignment information including columns:
#'   `qname`, `flag`, `rname`, `pos`, `mapq`, `cigar`.
#'
#' @examples
#' \donttest{
#' # Requires a BAM file
#' # bam_df <- bb_read_bam("aligned.bam")
#' }
#'
#' @export
bb_read_bam <- function(path, index = TRUE,
                        what = c("qname", "flag", "rname", "pos",
                                 "mapq", "cigar")) {

  if (!file.exists(path)) {
    stop("File not found: ", path, call. = FALSE)
  }

  if (requireNamespace("Rsamtools", quietly = TRUE)) {
    .read_bam_rsamtools(path, index, what)
  } else {
    message("Rsamtools not available. Trying system samtools.")
    supported <- c("qname", "flag", "rname", "pos", "mapq", "cigar")
    if (any(!what %in% supported)) {
      stop("System BAM reader supports: ", paste(supported, collapse = ", "), call. = FALSE)
    }
    .read_bam_system(path)[, what, drop = FALSE]
  }
}


#' Count Reads in BAM File
#'
#' Returns the total number of reads in a BAM file using Rsamtools or
#' system samtools.
#'
#' @param path Character. Path to a BAM file.
#'
#' @return A numeric value: total number of reads.
#'
#' @examples
#' \donttest{
#' # count <- bb_count_bam("aligned.bam")
#' }
#'
#' @export
bb_count_bam <- function(path) {

  if (!file.exists(path)) {
    stop("File not found: ", path, call. = FALSE)
  }

  if (requireNamespace("Rsamtools", quietly = TRUE)) {
    cnt <- Rsamtools::countBam(path)
    cnt$records
  } else {
    check_tool("samtools")
    out <- .run_tool("samtools", c("view", "-c", path))
    if (length(out) != 1L || !grepl("^[0-9]+$", trimws(out))) {
      stop("Invalid read count returned by samtools.", call. = FALSE)
    }
    as.numeric(trimws(out))
  }
}


#' @noRd
.read_bam_rsamtools <- function(path, index, what) {
  if (index) {
    bai <- paste0(path, ".bai")
    if (!file.exists(bai)) {
      Rsamtools::indexBam(path)
    }
  }
  param <- Rsamtools::ScanBamParam(what = what)
  bam <- Rsamtools::scanBam(path, param = param)[[1]]

  # Convert list to data.frame
  lens <- vapply(bam, length, integer(1))
  if (all(lens == 0L)) {
    return(data.frame(stringsAsFactors = FALSE))
  }

  as.data.frame(bam, stringsAsFactors = FALSE)
}


#' @noRd
.read_bam_system <- function(path) {
  check_tool("samtools")

  out <- .run_tool("samtools", c("view", path))

  if (length(out) == 0L) {
    return(data.frame(
      qname = character(0), flag = integer(0), rname = character(0),
      pos = integer(0), mapq = integer(0), cigar = character(0),
      stringsAsFactors = FALSE
    ))
  }

  # Parse SAM text output (tab-delimited)
  fields <- strsplit(out, "\t", fixed = TRUE)
  if (any(lengths(fields) < 11L)) {
    stop("Invalid SAM record returned by samtools.", call. = FALSE)
  }

  data.frame(
    qname = vapply(fields, `[`, character(1), 1L),
    flag  = as.integer(vapply(fields, `[`, character(1), 2L)),
    rname = vapply(fields, `[`, character(1), 3L),
    pos   = as.integer(vapply(fields, `[`, character(1), 4L)),
    mapq  = as.integer(vapply(fields, `[`, character(1), 5L)),
    cigar = vapply(fields, `[`, character(1), 6L),
    stringsAsFactors = FALSE
  )
}
