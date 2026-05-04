#' Generate QC Metrics
#'
#' Computes quality control metrics from FASTQ and/or BAM files.
#'
#' @param fastq_path Character or NULL. Path(s) to FASTQ file(s).
#' @param bam_path Character or NULL. Path(s) to BAM file(s).
#' @param use_fastqc Logical. Try to use system FastQC if available.
#'   Default `TRUE`.
#'
#' @return A `bb_qc` object containing:
#' \describe{
#'   \item{read_counts}{Total reads per file}
#'   \item{quality_scores}{Per-position quality score summary}
#'   \item{gc_content}{GC content distribution}
#'   \item{read_lengths}{Read length distribution}
#'   \item{mapping_rate}{Mapping rate from BAM files (if provided)}
#' }
#'
#' @examples
#' \donttest{
#' # Create a test FASTQ
#' tmp <- tempfile(fileext = ".fastq")
#' writeLines(c(
#'   "@read1", "ACGTACGT", "+", "IIIIIIII",
#'   "@read2", "GCGCGCGC", "+", "HHHHHHHH"
#' ), tmp)
#' qc <- bb_qc(fastq_path = tmp)
#' qc
#' }
#'
#' @export
bb_qc <- function(fastq_path = NULL, bam_path = NULL, use_fastqc = TRUE) {

  if (is.null(fastq_path) && is.null(bam_path)) {
    stop("Provide at least one of 'fastq_path' or 'bam_path'.", call. = FALSE)
  }

  result <- list(
    read_counts    = list(),
    quality_scores = list(),
    gc_content     = list(),
    read_lengths   = list(),
    mapping_rate   = list()
  )

  # FASTQ metrics
  if (!is.null(fastq_path)) {
    for (fq in fastq_path) {
      if (!file.exists(fq)) {
        warning("FASTQ file not found, skipping: ", fq, call. = FALSE)
        next
      }
      reads <- bb_read_fastq(fq)
      fname <- basename(fq)

      result$read_counts[[fname]] <- nrow(reads)
      result$gc_content[[fname]] <- .compute_gc(reads$sequence)
      result$read_lengths[[fname]] <- nchar(reads$sequence)
      result$quality_scores[[fname]] <- .compute_qual_summary(reads$quality)
    }
  }

  # BAM metrics
  if (!is.null(bam_path)) {
    for (bam in bam_path) {
      if (!file.exists(bam)) {
        warning("BAM file not found, skipping: ", bam, call. = FALSE)
        next
      }
      fname <- basename(bam)
      result$read_counts[[fname]] <- bb_count_bam(bam)
      result$mapping_rate[[fname]] <- .compute_mapping_rate(bam)
    }
  }

  structure(result, class = "bb_qc")
}


#' Print method for bb_qc
#' @param x A bb_qc object.
#' @param ... Additional arguments (ignored).
#' @export
print.bb_qc <- function(x, ...) {
  cat("bambamR QC Summary\n")
  cat("==================\n")
  cat("Files analyzed:", length(x$read_counts), "\n")
  for (nm in names(x$read_counts)) {
    cat(sprintf("  %s: %s reads", nm,
                format(x$read_counts[[nm]], big.mark = ",")))
    if (nm %in% names(x$mapping_rate)) {
      cat(sprintf(" (%.1f%% mapped)", x$mapping_rate[[nm]] * 100))
    }
    cat("\n")
  }
  invisible(x)
}


#' QC Summary
#'
#' Returns a summary data.frame of QC metrics.
#'
#' @param qc_object A `bb_qc` object from [bb_qc()].
#'
#' @return A data.frame with one row per file.
#'
#' @examples
#' \donttest{
#' tmp <- tempfile(fileext = ".fastq")
#' writeLines(c("@r1", "ACGT", "+", "IIII"), tmp)
#' qc <- bb_qc(fastq_path = tmp)
#' bb_qc_summary(qc)
#' }
#'
#' @export
bb_qc_summary <- function(qc_object) {
  if (!methods::is(qc_object, "bb_qc")) {
    stop("'qc_object' must be a bb_qc object.", call. = FALSE)
  }

  files <- names(qc_object$read_counts)
  data.frame(
    file = files,
    total_reads = vapply(files, function(f) qc_object$read_counts[[f]],
                         numeric(1)),
    median_gc = vapply(files, function(f) {
      gc <- qc_object$gc_content[[f]]
      if (is.null(gc)) NA_real_ else stats::median(gc)
    }, numeric(1)),
    mapping_rate = vapply(files, function(f) {
      mr <- qc_object$mapping_rate[[f]]
      if (is.null(mr)) NA_real_ else mr
    }, numeric(1)),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

#' Summary method for bb_qc
#' @param object A bb_qc object.
#' @param ... Additional arguments (ignored).
#' @export
summary.bb_qc <- function(object, ...) {
  bb_qc_summary(object)
}


# Internal QC helpers
# ==================================================================

#' Compute GC content per read
#' @noRd
.compute_gc <- function(sequences) {
  vapply(sequences, function(s) {
    chars <- strsplit(toupper(s), "", fixed = TRUE)[[1]]
    sum(chars %in% c("G", "C")) / length(chars)
  }, numeric(1), USE.NAMES = FALSE)
}


#' Compute per-position quality summary
#' @noRd
.compute_qual_summary <- function(quality_strings) {
  if (length(quality_strings) == 0L) return(NULL)

  # Convert Phred+33 to numeric
  qual_nums <- lapply(quality_strings, function(q) {
    as.integer(charToRaw(q)) - 33L
  })

  max_len <- max(vapply(qual_nums, length, integer(1)))

  # Per-position mean quality
  pos_means <- vapply(seq_len(max_len), function(pos) {
    vals <- vapply(qual_nums, function(q) {
      if (pos <= length(q)) q[pos] else NA_integer_
    }, integer(1))
    mean(vals, na.rm = TRUE)
  }, numeric(1))

  data.frame(
    position = seq_len(max_len),
    mean_quality = pos_means,
    stringsAsFactors = FALSE
  )
}


#' Compute mapping rate from BAM flagstat
#' @noRd
.compute_mapping_rate <- function(bam_path) {
  if (requireNamespace("Rsamtools", quietly = TRUE)) {
    cnt <- Rsamtools::countBam(bam_path)
    total <- cnt$records
    if (total == 0L) return(0)
    # Mapped reads: use flag filter
    param <- Rsamtools::ScanBamParam(
      flag = Rsamtools::scanBamFlag(isUnmappedQuery = FALSE)
    )
    mapped <- Rsamtools::countBam(bam_path, param = param)$records
    return(mapped / total)
  }

  # Fallback: system samtools
  if (nchar(Sys.which("samtools")) > 0L) {
    out <- system2("samtools", args = c("flagstat", bam_path),
                   stdout = TRUE, stderr = FALSE)
    # Parse first line: "N + N mapped"
    mapped_line <- grep("mapped", out, value = TRUE)[1]
    if (is.na(mapped_line)) return(NA_real_)
    total_line <- out[1]
    total <- as.numeric(sub(" \\+.*", "", total_line))
    mapped <- as.numeric(sub(" \\+.*", "", mapped_line))
    if (total == 0) return(0)
    return(mapped / total)
  }

  NA_real_
}
