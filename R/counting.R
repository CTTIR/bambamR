#' Count Reads per Gene/Feature
#'
#' Counts reads overlapping genomic features. Uses
#' `GenomicAlignments::summarizeOverlaps()` if available, otherwise shells out
#' to `featureCounts` from the Subread package.
#'
#' @param bam_paths Character vector. Paths to BAM files.
#' @param annotation Character. Path to GTF/GFF annotation file.
#' @param method Character. `"auto"` tries GenomicAlignments first, then
#'   featureCounts; `"internal"` forces GenomicAlignments; `"featureCounts"`
#'   forces the external tool.
#' @param threads Integer. Number of threads for featureCounts. Default `4`.
#' @param feature_type Character. Feature type to count (GTF column 3).
#'   Default `"exon"`.
#' @param attr_type Character. Attribute to group by. Default `"gene_id"`.
#' @param paired Logical. Paired-end data. Default `FALSE`.
#'
#' @return A numeric matrix of counts (genes x samples). Rownames are gene
#'   IDs, colnames are derived from BAM file names.
#'
#' @examples
#' \donttest{
#' # Requires BAM files and a GTF annotation
#' # counts <- bb_count_reads(c("s1.bam", "s2.bam"), "genes.gtf")
#' }
#'
#' @export
bb_count_reads <- function(bam_paths, annotation,
                           method = c("auto", "internal", "featureCounts"),
                           threads = 4L,
                           feature_type = "exon",
                           attr_type = "gene_id",
                           paired = FALSE) {

  method <- match.arg(method)

  # Validate BAM paths
  missing_bams <- bam_paths[!file.exists(bam_paths)]
  if (length(missing_bams) > 0L) {
    stop("BAM file(s) not found: ", paste(missing_bams, collapse = ", "),
         call. = FALSE)
  }
  if (!file.exists(annotation)) {
    stop("Annotation file not found: ", annotation, call. = FALSE)
  }

  if (method == "auto") {
    if (requireNamespace("GenomicAlignments", quietly = TRUE) &&
        requireNamespace("GenomicRanges", quietly = TRUE)) {
      return(.count_bioc(bam_paths, annotation, feature_type, attr_type,
                         paired))
    }
    if (nchar(Sys.which("featureCounts")) > 0L) {
      return(.count_featurecounts(bam_paths, annotation, threads,
                                  feature_type, attr_type, paired))
    }
    stop(
      "No counting method available. Install GenomicAlignments or featureCounts.",
      call. = FALSE
    )
  }

  if (method == "internal") {
    check_pkg("GenomicAlignments", reason = "for read counting")
    check_pkg("GenomicRanges", reason = "for read counting")
    return(.count_bioc(bam_paths, annotation, feature_type, attr_type, paired))
  }

  if (method == "featureCounts") {
    check_tool("featureCounts")
    return(.count_featurecounts(bam_paths, annotation, threads,
                                feature_type, attr_type, paired))
  }
}


# Bioconductor-based counting
# ------------------------------------------------------------------
#' @noRd
.count_bioc <- function(bam_paths, annotation, feature_type, attr_type,
                        paired) {

  check_pkg("GenomicFeatures", reason = "for GTF parsing")
  check_pkg("Rsamtools", reason = "for BAM reading")

  message("Counting reads with GenomicAlignments...")

  # Load annotation
  txdb <- GenomicFeatures::makeTxDbFromGFF(annotation, format = "gtf")
  features <- GenomicFeatures::exonsBy(txdb, by = "gene")

  # BAM files
  bfl <- Rsamtools::BamFileList(bam_paths)

  # Count
  mode <- if (paired) "IntersectionNotEmpty" else "Union"
  se <- GenomicAlignments::summarizeOverlaps(
    features = features,
    reads = bfl,
    mode = mode,
    singleEnd = !paired,
    ignore.strand = TRUE
  )

  counts <- SummarizedExperiment::assay(se)
  colnames(counts) <- .clean_sample_names(bam_paths)
  counts
}


# featureCounts-based counting
# ------------------------------------------------------------------
#' @noRd
.count_featurecounts <- function(bam_paths, annotation, threads,
                                 feature_type, attr_type, paired) {

  check_tool("featureCounts")
  message("Counting reads with featureCounts...")

  out_file <- tempfile(fileext = ".txt")

  args <- c(
    "-a", annotation,
    "-o", out_file,
    "-t", feature_type,
    "-g", attr_type,
    "-T", threads
  )
  if (paired) args <- c(args, "-p")
  args <- c(args, bam_paths)

  system2("featureCounts", args = args, stdout = FALSE, stderr = FALSE)

  if (!file.exists(out_file)) {
    stop("featureCounts failed to produce output.", call. = FALSE)
  }

  # Parse output
  result <- data.table::fread(out_file, skip = 1, header = TRUE,
                               data.table = FALSE)
  gene_ids <- result[[1]]  # Geneid column
  count_cols <- seq(7, ncol(result))  # Count columns start at column 7
  counts <- as.matrix(result[, count_cols, drop = FALSE])
  rownames(counts) <- gene_ids
  colnames(counts) <- .clean_sample_names(bam_paths)

  unlink(out_file)
  unlink(paste0(out_file, ".summary"))

  counts
}


#' Clean sample names from BAM paths
#' @noRd
.clean_sample_names <- function(bam_paths) {
  nms <- basename(bam_paths)
  sub("\\.[Bb][Aa][Mm]$", "", nms)
}
