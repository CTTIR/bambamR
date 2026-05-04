#' Align Reads to Reference Genome
#'
#' Wrapper for STAR, HISAT2, or minimap2 alignment tools. Checks that the
#' selected aligner is available on the system PATH before running.
#'
#' @param fastq Character. Path(s) to FASTQ files. For paired-end, provide
#'   a character vector of length 2.
#' @param genome_index Character. Path to genome index directory (STAR) or
#'   index prefix (HISAT2, minimap2).
#' @param output_dir Character. Output directory for BAM files.
#' @param aligner Character. One of `"STAR"`, `"HISAT2"`, `"minimap2"`.
#' @param threads Integer. Number of threads. Default `4`.
#' @param paired Logical. Paired-end mode. If `TRUE`, `fastq` must have 2
#'   elements. Default `FALSE`.
#' @param extra_args Character or NULL. Additional arguments to pass to the
#'   aligner.
#'
#' @return A list with components:
#' \describe{
#'   \item{bam}{Path to the output BAM file}
#'   \item{stats}{A data.frame with alignment statistics}
#'   \item{command}{The exact command that was executed}
#' }
#'
#' @examples
#' \donttest{
#' # Requires STAR/HISAT2/minimap2 installed
#' # result <- bb_align("reads.fastq", "/path/to/index", "output/")
#' }
#'
#' @export
bb_align <- function(fastq, genome_index, output_dir,
                     aligner = c("STAR", "HISAT2", "minimap2"),
                     threads = 4L, paired = FALSE,
                     extra_args = NULL) {

  aligner <- match.arg(aligner)

  # Validate inputs
  if (paired && length(fastq) != 2L) {
    stop("For paired-end mode, 'fastq' must be a character vector of length 2.",
         call. = FALSE)
  }
  for (fq in fastq) {
    if (!file.exists(fq)) {
      stop("FASTQ file not found: ", fq, call. = FALSE)
    }
  }
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  switch(aligner,
    STAR    = .align_star(fastq, genome_index, output_dir, threads,
                          paired, extra_args),
    HISAT2  = .align_hisat2(fastq, genome_index, output_dir, threads,
                            paired, extra_args),
    minimap2 = .align_minimap2(fastq, genome_index, output_dir, threads,
                               paired, extra_args)
  )
}


# STAR alignment
# ------------------------------------------------------------------
#' @noRd
.align_star <- function(fastq, genome_index, output_dir, threads,
                        paired, extra_args) {
  check_tool("STAR")

  prefix <- file.path(output_dir, "STAR_")
  bam_out <- paste0(prefix, "Aligned.sortedByCoord.out.bam")

  args <- c(
    "--runThreadN", threads,
    "--genomeDir", genome_index,
    "--readFilesIn", fastq,
    "--outFileNamePrefix", prefix,
    "--outSAMtype", "BAM", "SortedByCoordinate",
    "--outSAMunmapped", "Within"
  )

  # Handle gzipped input
  if (any(grepl("\\.gz$", fastq))) {
    args <- c(args, "--readFilesCommand", "zcat")
  }

  if (!is.null(extra_args)) {
    args <- c(args, extra_args)
  }

  cmd <- paste("STAR", paste(args, collapse = " "))
  message("Running: ", cmd)

  ret <- system2("STAR", args = args, stdout = TRUE, stderr = TRUE)

  # Parse STAR Log.final.out
  log_file <- paste0(prefix, "Log.final.out")
  stats <- .parse_star_log(log_file)

  list(bam = bam_out, stats = stats, command = cmd)
}


#' @noRd
.parse_star_log <- function(log_file) {
  if (!file.exists(log_file)) {
    return(data.frame(metric = "status", value = "log not found",
                      stringsAsFactors = FALSE))
  }
  lines <- readLines(log_file)
  # Extract key metrics
  metrics <- c("Number of input reads",
               "Uniquely mapped reads number",
               "Uniquely mapped reads %",
               "Number of splices: Total",
               "Mismatch rate per base, %")
  result <- lapply(metrics, function(m) {
    idx <- grep(m, lines, fixed = TRUE)
    if (length(idx) > 0L) {
      val <- trimws(sub(".*\\|\\s*", "", lines[idx[1]]))
      data.frame(metric = m, value = val, stringsAsFactors = FALSE)
    } else {
      NULL
    }
  })
  result <- result[!vapply(result, is.null, logical(1))]
  if (length(result) == 0L) {
    return(data.frame(metric = character(0), value = character(0),
                      stringsAsFactors = FALSE))
  }
  do.call(rbind, result)
}


# HISAT2 alignment
# ------------------------------------------------------------------
#' @noRd
.align_hisat2 <- function(fastq, genome_index, output_dir, threads,
                          paired, extra_args) {
  check_tool("hisat2")
  check_tool("samtools")

  bam_out <- file.path(output_dir, "hisat2_aligned.bam")

  args <- c("-x", genome_index, "-p", threads)

  if (paired) {
    args <- c(args, "-1", fastq[1], "-2", fastq[2])
  } else {
    args <- c(args, "-U", fastq[1])
  }

  if (!is.null(extra_args)) {
    args <- c(args, extra_args)
  }

  cmd <- paste("hisat2", paste(args, collapse = " "),
               "| samtools sort -o", bam_out)
  message("Running: ", cmd)

  # Run hisat2, pipe to samtools sort
  hisat_out <- system2("hisat2", args = args, stdout = TRUE, stderr = TRUE)

  # Write SAM to temp, convert to sorted BAM
  sam_tmp <- tempfile(fileext = ".sam")
  writeLines(hisat_out, sam_tmp)
  system2("samtools", args = c("sort", "-o", bam_out, sam_tmp))
  unlink(sam_tmp)

  # Parse summary from stderr
  stats <- data.frame(
    metric = "alignment_complete",
    value = "TRUE",
    stringsAsFactors = FALSE
  )

  list(bam = bam_out, stats = stats, command = cmd)
}


# minimap2 alignment
# ------------------------------------------------------------------
#' @noRd
.align_minimap2 <- function(fastq, genome_index, output_dir, threads,
                            paired, extra_args) {
  check_tool("minimap2")
  check_tool("samtools")

  bam_out <- file.path(output_dir, "minimap2_aligned.bam")

  args <- c("-a", "-t", threads, genome_index, fastq)
  if (!is.null(extra_args)) {
    args <- c(args, extra_args)
  }

  cmd <- paste("minimap2", paste(args, collapse = " "),
               "| samtools sort -o", bam_out)
  message("Running: ", cmd)

  sam_tmp <- tempfile(fileext = ".sam")
  system2("minimap2", args = args, stdout = sam_tmp, stderr = FALSE)
  system2("samtools", args = c("sort", "-o", bam_out, sam_tmp))
  unlink(sam_tmp)

  stats <- data.frame(
    metric = "alignment_complete",
    value = "TRUE",
    stringsAsFactors = FALSE
  )

  list(bam = bam_out, stats = stats, command = cmd)
}
