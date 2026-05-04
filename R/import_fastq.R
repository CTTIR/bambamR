#' Read FASTQ Files
#'
#' Reads FASTQ files using ShortRead (if available) or a minimal
#' base-R parser as fallback. Supports `.gz` compressed files.
#'
#' @param path Character. Path to a FASTQ or FASTQ.gz file.
#' @param n Integer or NULL. Number of reads to sample. `NULL` reads all.
#'
#' @return A data.frame with columns:
#' \describe{
#'   \item{id}{Read identifier}
#'   \item{sequence}{Nucleotide sequence}
#'   \item{quality}{Quality string (Phred+33 encoded)}
#' }
#'
#' @examples
#' \donttest{
#' # Create a temporary FASTQ file
#' tmp <- tempfile(fileext = ".fastq")
#' writeLines(c(
#'   "@read1", "ACGTACGT", "+", "IIIIIIII",
#'   "@read2", "TGCATGCA", "+", "HHHHHHHH"
#' ), tmp)
#' reads <- bb_read_fastq(tmp)
#' reads
#' }
#'
#' @export
bb_read_fastq <- function(path, n = NULL) {

  if (!file.exists(path)) {
    stop("File not found: ", path, call. = FALSE)
  }

  if (requireNamespace("ShortRead", quietly = TRUE)) {
    .read_fastq_shortread(path, n)
  } else {
    message("ShortRead not available. Using base-R FASTQ parser.")
    .read_fastq_base(path, n)
  }
}


#' Read FASTQ with ShortRead
#' @noRd
.read_fastq_shortread <- function(path, n) {
  fq <- ShortRead::readFastq(path)
  if (!is.null(n)) {
    n <- min(n, length(fq))
    fq <- fq[seq_len(n)]
  }
  # quality() generic lives in Biostrings; the returned FastqQuality object
  # stores the raw strings in an internal "quality" slot (a BStringSet).
  qual_obj <- Biostrings::quality(fq)
  qual_strings <- as.character(methods::slot(qual_obj, "quality"))
  data.frame(
    id       = as.character(ShortRead::id(fq)),
    sequence = as.character(ShortRead::sread(fq)),
    quality  = qual_strings,
    stringsAsFactors = FALSE
  )
}


#' Read FASTQ with base R (fallback)
#' @noRd
.read_fastq_base <- function(path, n) {
  # Handle gzipped files
  if (grepl("\\.gz$", path)) {
    con <- gzfile(path, "r")
    on.exit(close(con))
    lines <- readLines(con)
  } else {
    lines <- readLines(path)
  }

  # FASTQ has 4 lines per record
  n_records <- length(lines) %/% 4L
  if (n_records == 0L) {
    return(data.frame(id = character(0), sequence = character(0),
                      quality = character(0), stringsAsFactors = FALSE))
  }

  if (!is.null(n)) {
    n_records <- min(n, n_records)
  }

  idx <- seq_len(n_records)
  line1 <- (idx - 1L) * 4L + 1L  # ID lines
  line2 <- line1 + 1L             # Sequence lines
  line4 <- line1 + 3L             # Quality lines

  data.frame(
    id       = sub("^@", "", lines[line1]),
    sequence = lines[line2],
    quality  = lines[line4],
    stringsAsFactors = FALSE
  )
}
