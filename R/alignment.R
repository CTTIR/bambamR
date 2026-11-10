#' Align Reads to Reference Genome
#'
#' Wrapper for STAR, HISAT2, or minimap2 alignment tools. Checks that the
#' selected aligner is available on the system PATH before running.
#'
#' @param fastq Character. Path(s) to FASTQ files. For paired-end, provide
#'   a character vector of length 2.
#' @param genome_index Character. Path to genome index directory (STAR) or
#'   index prefix (HISAT2, minimap2).
#' @param output_dir Character. New or empty output directory for BAM files.
#' @param aligner Character. One of `"STAR"`, `"HISAT2"`, `"minimap2"`.
#' @param threads Integer. Number of threads. Default `4`.
#' @param paired Logical. Paired-end mode. If `TRUE`, `fastq` must have 2
#'   elements. Default `FALSE`.
#' @param extra_args Character or NULL. Additional arguments to pass to the
#'   aligner, one argument per element (for example `c("-x", "sr")`).
#'   HISAT2 optional arguments must not contain whitespace or shell metacharacters:
#'   its upstream wrapper reparses them through a shell. FASTQ and index paths
#'   containing these characters are handled with temporary symbolic links.
#'
#' @return A list with components:
#' \describe{
#'   \item{bam}{Path to the output BAM file}
#'   \item{stats}{A data.frame with alignment statistics}
#'   \item{command}{The exact command that was executed}
#'   \item{input_paths}{For HISAT2, the original FASTQ and index paths before
#'     temporary aliasing. Temporary aliases are removed after execution.}
#' }
#'
#' @examples
#' \dontrun{
#' # Not run: needs an external aligner (STAR, HISAT2 or minimap2) on the
#' # PATH and a prebuilt genome index, neither of which can be shipped.
#' result <- bb_align("reads.fastq", "/path/to/index",
#'                    file.path(tempdir(), "alignment"), aligner = "HISAT2")
#' }
#'
#' @export
bb_align <- function(fastq, genome_index, output_dir,
                     aligner = c("STAR", "HISAT2", "minimap2"),
                     threads = 4L, paired = FALSE,
                     extra_args = NULL) {

  aligner <- match.arg(aligner)
  .validate_threads(threads)
  if (!is.logical(paired) || length(paired) != 1L || is.na(paired)) {
    stop("'paired' must be TRUE or FALSE.", call. = FALSE)
  }

  # Validate inputs
  if (paired && length(fastq) != 2L) {
    stop("For paired-end mode, 'fastq' must be a character vector of length 2.",
         call. = FALSE)
  }
  if (!is.character(fastq) || anyNA(fastq) || (!paired && length(fastq) != 1L)) {
    stop("Single-end mode requires one FASTQ path.", call. = FALSE)
  }
  for (fq in fastq) {
    if (!file.exists(fq)) {
      stop("FASTQ file not found: ", fq, call. = FALSE)
    }
  }
  if (!dir.exists(output_dir)) {
    if (!dir.create(output_dir, recursive = TRUE)) {
      stop("Cannot create output directory: ", output_dir, call. = FALSE)
    }
  } else if (length(list.files(output_dir, all.files = TRUE, no.. = TRUE))) {
    stop("Alignment output directory must be empty: ", output_dir, call. = FALSE)
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

  cmd <- paste("STAR", paste(shQuote(as.character(args)), collapse = " "))
  message("Running: ", cmd)

  .run_tool("STAR", args)
  .require_product(bam_out)

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

  if (!is.null(extra_args) && (!is.character(extra_args) || anyNA(extra_args) ||
      any(!grepl("^[A-Za-z0-9_./,:=+%-]+$", extra_args)))) {
    stop("HISAT2 optional arguments cannot safely contain whitespace or shell metacharacters; ",
         "its upstream wrapper reparses them through a shell.", call. = FALSE)
  }
  input_paths <- list(fastq = fastq, genome_index = genome_index)
  inputs <- .hisat2_inputs(fastq, genome_index)
  if (!is.null(inputs$directory)) on.exit(unlink(inputs$directory, recursive = TRUE), add = TRUE)
  fastq <- inputs$fastq
  genome_index <- inputs$genome_index

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

  sam_tmp <- tempfile(fileext = ".sam")
  on.exit(unlink(sam_tmp), add = TRUE)
  cmd <- paste("hisat2", paste(shQuote(as.character(args)), collapse = " "),
               ">", shQuote(sam_tmp), "&& samtools sort -o", shQuote(bam_out), shQuote(sam_tmp))
  message("Running: ", cmd)
  .run_tool("hisat2", args, stdout = sam_tmp)
  .require_product(sam_tmp)
  .run_tool("samtools", c("sort", "-o", bam_out, sam_tmp))
  .require_product(bam_out)

  # Parse summary from stderr
  stats <- data.frame(
    metric = "alignment_complete",
    value = "TRUE",
    stringsAsFactors = FALSE
  )

  list(bam = bam_out, stats = stats, command = cmd, input_paths = input_paths)
}

# The HISAT2 2.2.1 Perl wrapper joins arguments into another shell command.
# Quoting at system2's boundary alone therefore cannot preserve these paths.
# Alias only inputs, never copy large indexes or rewrite the user's files.
.hisat2_inputs <- function(fastq, genome_index) {
  # Backslashes are Windows path separators, but literal filename characters on Unix.
  if (.Platform$OS.type == "windows") {
    fastq <- gsub("\\", "/", fastq, fixed = TRUE)
    existing <- file.exists(fastq)
    fastq[existing] <- normalizePath(fastq[existing], winslash = "/", mustWork = TRUE)
    genome_index <- gsub("\\", "/", genome_index, fixed = TRUE)
  }
  safe <- function(x) grepl("^[A-Za-z0-9_./:+-]+$", x)
  if (all(safe(c(fastq, genome_index)))) {
    return(list(fastq = fastq, genome_index = genome_index, directory = NULL))
  }
  temporary_root <- normalizePath(tempdir(), winslash = "/", mustWork = TRUE)
  if (!safe(temporary_root) && .Platform$OS.type == "unix") temporary_root <- "/tmp"
  if (!safe(temporary_root)) stop("HISAT2 requires a temporary directory without shell metacharacters.",
                                  call. = FALSE)
  directory <- tempfile("bambamR-hisat2-", tmpdir = temporary_root)
  if (!dir.create(directory, mode = "0700")) stop("Cannot create HISAT2 staging directory.", call. = FALSE)
  complete <- FALSE
  on.exit(if (!complete) unlink(directory, recursive = TRUE))
  prefixes <- genome_index
  search_path <- Sys.getenv("HISAT2_INDEXES")
  if (nzchar(search_path)) prefixes <- c(prefixes, file.path(search_path, genome_index))
  suffixes <- character()
  for (prefix in prefixes) {
    for (extension in c("ht2", "ht2l")) {
      parts <- paste0(".", seq_len(8L), ".", extension)
      if (all(file.exists(paste0(prefix, parts)))) suffixes <- c(suffixes, parts)
    }
    if (length(suffixes)) break
  }
  if (!length(suffixes)) stop("Cannot find a complete HISAT2 index for: ", genome_index, call. = FALSE)
  sources <- c(paste0(prefix, suffixes), fastq)
  compression <- ifelse(grepl("\\.gz$", fastq), ".gz",
                         ifelse(grepl("\\.bz2$", fastq), ".bz2", ""))
  aliases <- file.path(directory, c(paste0("index", suffixes),
    paste0("reads", seq_along(fastq), ".fastq", compression)))
  sources <- normalizePath(sources, winslash = "/", mustWork = TRUE)
  if (!all(file.symlink(sources, aliases))) stop("Cannot create safe HISAT2 input aliases.", call. = FALSE)
  complete <- TRUE
  list(fastq = utils::tail(aliases, length(fastq)),
       genome_index = file.path(directory, "index"), directory = directory)
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

  sam_tmp <- tempfile(fileext = ".sam")
  on.exit(unlink(sam_tmp), add = TRUE)
  cmd <- paste("minimap2", paste(shQuote(as.character(args)), collapse = " "),
               ">", shQuote(sam_tmp), "&& samtools sort -o", shQuote(bam_out), shQuote(sam_tmp))
  message("Running: ", cmd)
  .run_tool("minimap2", args, stdout = sam_tmp)
  .require_product(sam_tmp)
  .run_tool("samtools", c("sort", "-o", bam_out, sam_tmp))
  .require_product(bam_out)

  stats <- data.frame(
    metric = "alignment_complete",
    value = "TRUE",
    stringsAsFactors = FALSE
  )

  list(bam = bam_out, stats = stats, command = cmd)
}
