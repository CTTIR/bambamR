#' Run Full bambamR Pipeline
#'
#' Orchestrates the complete RNA-seq analysis pipeline from FASTQ/BAM files
#' or a count matrix through alignment, counting, normalization, differential
#' expression, and visualization.
#'
#' @param fastq_dir Character or NULL. Directory containing FASTQ files.
#'   If provided, the pipeline starts from alignment.
#' @param bam_dir Character or NULL. Directory containing BAM files.
#'   If provided, the pipeline starts from read counting.
#' @param count_matrix A numeric matrix or NULL. Pre-computed count matrix.
#'   If provided, the pipeline starts from normalization/DE.
#' @param output_dir Character. Output directory for results.
#'   Default `"bambamR_output"`.
#' @param genome_index Character or NULL. Path to genome index. Required if
#'   starting from FASTQ.
#' @param annotation Character or NULL. Path to GTF annotation. Required if
#'   starting from FASTQ or BAM.
#' @param sample_info A data.frame with sample metadata. Must contain a
#'   `condition` column (or the variable in `design`). Rownames should
#'   match sample identifiers.
#' @param aligner Character. Aligner for FASTQ alignment. Default `"STAR"`.
#' @param de_method Character. DE method: `"DESeq2"`, `"edgeR"`, or
#'   `"limma"`. Default `"DESeq2"`.
#' @param design A formula for DE analysis. Default `~ condition`.
#' @param skip Character vector. Steps to skip: `"qc"`, `"align"`,
#'   `"count"`, `"de"`, `"viz"`. Default `character(0)`.
#' @param threads Integer. Number of threads. Default `4`.
#'
#' @return A `bb_result` object containing counts, metadata, DE results,
#'   and plots.
#'
#' @details
#' The pipeline allows entry at any stage:
#' - **From FASTQ**: requires `fastq_dir`, `genome_index`, `annotation`,
#'   and `sample_info`
#' - **From BAM**: requires `bam_dir`, `annotation`, and `sample_info`
#' - **From counts**: requires `count_matrix` and `sample_info`
#'
#' @examples
#' \donttest{
#' # Starting from a count matrix (minimal mode, no Bioconductor needed)
#' set.seed(42)
#' counts <- matrix(rpois(600, 100), nrow = 100, ncol = 6,
#'   dimnames = list(paste0("gene", 1:100), paste0("S", 1:6)))
#' sample_info <- data.frame(
#'   condition = factor(rep(c("ctrl", "treat"), each = 3)),
#'   row.names = paste0("S", 1:6)
#' )
#' # result <- bb_pipeline(count_matrix = counts, sample_info = sample_info)
#' }
#'
#' @export
bb_pipeline <- function(fastq_dir = NULL,
                        bam_dir = NULL,
                        count_matrix = NULL,
                        output_dir = "bambamR_output",
                        genome_index = NULL,
                        annotation = NULL,
                        sample_info = NULL,
                        aligner = "STAR",
                        de_method = c("DESeq2", "edgeR", "limma"),
                        design = ~ condition,
                        skip = character(0),
                        threads = 4L) {

  de_method <- match.arg(de_method)

  # Determine starting point
  start <- .determine_start(fastq_dir, bam_dir, count_matrix)

  # Validate inputs for the starting point
  .validate_pipeline_inputs(start, fastq_dir, bam_dir, count_matrix,
                            genome_index, annotation, sample_info)

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  result <- new_bb_result()
  result$metadata <- sample_info

  # Step 1: QC (FASTQ)
  if (start == "fastq" && !"qc" %in% skip) {
    message("=== Step 1: Quality Control ===")
    fq_files <- .list_fastq(fastq_dir)
    result$qc <- bb_qc(fastq_path = fq_files)
    message("  QC complete for ", length(fq_files), " files")
  }

  # Step 2: Alignment
  if (start == "fastq" && !"align" %in% skip) {
    message("=== Step 2: Alignment ===")
    fq_files <- .list_fastq(fastq_dir)
    align_dir <- file.path(output_dir, "alignment")
    align_results <- lapply(fq_files, function(fq) {
      bb_align(fq, genome_index, align_dir, aligner = aligner,
               threads = threads)
    })
    bam_dir <- align_dir
    result$alignment_stats <- do.call(rbind, lapply(align_results,
                                                     function(x) x$stats))
    message("  Alignment complete")
  }

  # Step 3: Read counting
  if (start %in% c("fastq", "bam") && !"count" %in% skip) {
    message("=== Step 3: Read Counting ===")
    bam_files <- list.files(bam_dir, pattern = "\\.bam$", full.names = TRUE)
    if (length(bam_files) == 0L) {
      stop("No BAM files found in: ", bam_dir, call. = FALSE)
    }
    count_matrix <- bb_count_reads(bam_files, annotation, threads = threads)
    message("  Counted ", nrow(count_matrix), " genes across ",
            ncol(count_matrix), " samples")
  }

  result$counts <- count_matrix

  # Step 4: Normalization
  message("=== Step 4: Normalization ===")
  norm_counts <- bb_normalize(count_matrix, method = "cpm")
  message("  CPM normalization complete")

  # Step 5: Differential Expression
  if (!"de" %in% skip && !is.null(sample_info)) {
    message("=== Step 5: Differential Expression (", de_method, ") ===")
    result$de_results <- tryCatch(
      .run_de(count_matrix, sample_info, de_method, design),
      error = function(e) {
        warning("DE analysis failed: ", conditionMessage(e), call. = FALSE)
        NULL
      }
    )
    if (!is.null(result$de_results)) {
      n_sig <- sum(result$de_results$padj < 0.05, na.rm = TRUE)
      message("  Found ", n_sig, " significant genes (FDR < 0.05)")
    }
  }

  # Step 6: Visualization
  if (!"viz" %in% skip) {
    message("=== Step 6: Generating Plots ===")
    result$plots <- list()

    # PCA
    if (!is.null(sample_info) && "condition" %in% colnames(sample_info)) {
      result$plots$pca <- tryCatch(
        bb_pca(norm_counts, sample_info, color_by = "condition"),
        error = function(e) NULL
      )
    }

    # Heatmap
    result$plots$heatmap <- tryCatch(
      bb_heatmap(norm_counts, de_result = result$de_results, n_genes = 50L),
      error = function(e) NULL
    )

    # Volcano + MA (only if DE results available)
    if (!is.null(result$de_results)) {
      result$plots$volcano <- tryCatch(
        bb_volcano(result$de_results),
        error = function(e) NULL
      )
      if ("basemean" %in% colnames(result$de_results)) {
        result$plots$ma <- tryCatch(
          bb_ma_plot(result$de_results),
          error = function(e) NULL
        )
      }
    }

    plot_names <- names(result$plots)[!vapply(result$plots, is.null,
                                               logical(1))]
    message("  Generated: ", paste(plot_names, collapse = ", "))
  }

  # Save intermediate results
  bb_export_rds(result, file.path(output_dir, "bambamR_result.rds"))

  message("\n=== Pipeline Complete ===")
  result
}


# Internal pipeline helpers
# ==================================================================

#' @noRd
.determine_start <- function(fastq_dir, bam_dir, count_matrix) {
  if (!is.null(fastq_dir)) return("fastq")
  if (!is.null(bam_dir)) return("bam")
  if (!is.null(count_matrix)) return("counts")
  stop("Provide one of: fastq_dir, bam_dir, or count_matrix.", call. = FALSE)
}

#' @noRd
.validate_pipeline_inputs <- function(start, fastq_dir, bam_dir,
                                      count_matrix, genome_index,
                                      annotation, sample_info) {
  if (start == "fastq") {
    if (!dir.exists(fastq_dir)) {
      stop("fastq_dir not found: ", fastq_dir, call. = FALSE)
    }
    if (is.null(genome_index)) {
      stop("'genome_index' is required when starting from FASTQ.", call. = FALSE)
    }
    if (is.null(annotation)) {
      stop("'annotation' is required when starting from FASTQ.", call. = FALSE)
    }
  }
  if (start == "bam") {
    if (!dir.exists(bam_dir)) {
      stop("bam_dir not found: ", bam_dir, call. = FALSE)
    }
    if (is.null(annotation)) {
      stop("'annotation' is required when starting from BAM.", call. = FALSE)
    }
  }
  if (start == "counts" && is.null(sample_info)) {
    stop("'sample_info' is required when starting from a count matrix.",
         call. = FALSE)
  }
}

#' @noRd
.list_fastq <- function(dir) {
  fq <- list.files(dir, pattern = "\\.(fastq|fq)(\\.gz)?$",
                   full.names = TRUE)
  if (length(fq) == 0L) {
    stop("No FASTQ files found in: ", dir, call. = FALSE)
  }
  fq
}

#' @noRd
.run_de <- function(counts, sample_info, method, design) {
  switch(method,
    DESeq2 = bb_deseq2(counts, sample_info, design = design),
    edgeR  = {
      if (!"condition" %in% colnames(sample_info)) {
        stop("'sample_info' must contain a 'condition' column for edgeR.",
             call. = FALSE)
      }
      shared <- intersect(colnames(counts), rownames(sample_info))
      bb_edger(counts[, shared, drop = FALSE],
               sample_info[shared, "condition"])
    },
    limma = {
      shared <- intersect(colnames(counts), rownames(sample_info))
      dm <- stats::model.matrix(design, data = sample_info[shared, ,
                                                            drop = FALSE])
      bb_limma_voom(counts[, shared, drop = FALSE], dm)
    }
  )
}
