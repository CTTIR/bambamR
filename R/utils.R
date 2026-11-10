# Internal utilities for bambamR
# ================================

#' @keywords internal
"_PACKAGE"

# Package-level global variables to avoid R CMD check NOTEs
utils::globalVariables(c(

  # data.table columns
  ".", ".N", ".SD",


  # oncoplot / viz columns
  "gene", "sample", "mutation_type", "n_mutations", "pct", "pct_label",
  "gene_count", "sample_count", "value", "annotation_value",


  # DE result columns
  "log2fc", "padj", "pvalue", "basemean", "sig", "label",

  # QC columns
  "position", "quality_score", "count", "gc_content", "read_length",
  "mean_quality", "file",

  # PCA columns
  "PC1", "PC2", "sample_name", ".data",

  # heatmap columns
  "variable", "row_gene", "col_sample"
))


# .onAttach hook (startup messages must go in .onAttach, not .onLoad)
# ------------------------------------------------------------------
.onAttach <- function(libname, pkgname) {

  bioc_pkgs <- c("DESeq2", "edgeR", "limma", "ShortRead", "Rsamtools",
                 "GenomicAlignments", "GenomicRanges", "SummarizedExperiment",
                 "GenomicFeatures", "txdbmaker", "ComplexHeatmap")
  available <- vapply(bioc_pkgs, function(pkg) nzchar(system.file(package = pkg)), logical(1))

  if (all(available)) {
    packageStartupMessage("bambamR: Optional Bioconductor engines are installed")
  } else {
    missing <- bioc_pkgs[!available]
    packageStartupMessage(
      "bambamR: Some optional engines are not installed. ",
      "Optional packages not found: ", paste(missing, collapse = ", "), "\n",
      "Install with: BiocManager::install(c(",
      paste0("'", missing, "'", collapse = ", "), "))"
    )
  }
}


# Internal helper: check if a package is available
# ------------------------------------------------------------------
#' Check if a package is available
#'
#' @param pkg Character. Package name.
#' @param reason Character. Why the package is needed (for error message).
#' @param install_cmd Character or NULL. Installation command hint.
#' @return Logical TRUE invisibly if available, otherwise stops with message.
#' @noRd
check_pkg <- function(pkg, reason = NULL, install_cmd = NULL) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    return(invisible(TRUE))
  }

  msg <- paste0("Package '", pkg, "' is required")
  if (!is.null(reason)) {
    msg <- paste0(msg, " ", reason)
  }
  msg <- paste0(msg, ".")

  if (is.null(install_cmd)) {
    # Guess installation method
    bioc_pkgs <- c("DESeq2", "edgeR", "limma", "ShortRead", "Rsamtools",
                   "GenomicAlignments", "GenomicRanges", "SummarizedExperiment",
                   "Biobase", "ComplexHeatmap", "GenomicFeatures", "txdbmaker")
    if (pkg %in% bioc_pkgs) {
      install_cmd <- paste0("BiocManager::install('", pkg, "')")
    } else {
      install_cmd <- paste0("install.packages('", pkg, "')")
    }
  }
  msg <- paste0(msg, "\nInstall with: ", install_cmd)
  stop(msg, call. = FALSE)
}


# Internal helper: check if a system tool is available
# ------------------------------------------------------------------
#' @noRd
check_tool <- function(tool) {
  path <- Sys.which(tool)
  if (nchar(path) == 0L) {
    stop("External tool '", tool, "' not found on PATH.", call. = FALSE)
  }
  invisible(path)
}


# Internal helper: validate count matrix
# ------------------------------------------------------------------
#' @noRd
validate_counts <- function(counts) {
  if (!is.matrix(counts) && !is.data.frame(counts)) {
    stop("'counts' must be a matrix or data.frame.", call. = FALSE)
  }
  if (is.data.frame(counts)) {
    counts <- as.matrix(counts)
  }
  if (!is.numeric(counts)) {
    stop("'counts' must contain numeric values.", call. = FALSE)
  }
  if (is.null(rownames(counts))) {
    stop("'counts' must have rownames (gene identifiers).", call. = FALSE)
  }
  if (is.null(colnames(counts))) {
    stop("'counts' must have colnames (sample identifiers).", call. = FALSE)
  }
  if (any(dim(counts) == 0L) || any(!is.finite(counts)) || any(counts < 0)) {
    stop("'counts' must be non-empty, finite and non-negative.", call. = FALSE)
  }
  for (ids in dimnames(counts)) {
    if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
      stop("Gene and sample identifiers must be non-empty and unique.", call. = FALSE)
    }
  }
  counts
}

#' @noRd
.match_metadata <- function(counts, metadata) {
  if (!is.data.frame(metadata) || is.null(rownames(metadata)) ||
        anyDuplicated(rownames(metadata)) ||
        !setequal(colnames(counts), rownames(metadata))) {
    stop("Metadata rownames must match all count sample names exactly.", call. = FALSE)
  }
  metadata[colnames(counts), , drop = FALSE]
}

#' @noRd
.validate_libraries <- function(counts) {
  totals <- colSums(counts)
  if (any(!is.finite(totals) | totals <= 0)) {
    stop("Every sample library must have a finite, positive total count.", call. = FALSE)
  }
  invisible(counts)
}

#' @noRd
.validate_design <- function(design, counts) {
  if (!is.matrix(design) || !is.numeric(design) ||
        nrow(design) != ncol(counts) || ncol(design) < 1L ||
        any(!is.finite(design))) {
    stop("'design' must be a finite numeric matrix with one row per sample.", call. = FALSE)
  }
  ids <- rownames(design)
  if (!is.null(ids) && !identical(ids, as.character(seq_len(nrow(design))))) {
    if (anyDuplicated(ids) || !setequal(ids, colnames(counts))) {
      stop("Design rownames must match count sample names.", call. = FALSE)
    }
    design <- design[colnames(counts), , drop = FALSE]
  }
  if (qr(design)$rank < ncol(design) || nrow(design) <= ncol(design)) {
    stop("Design must have full rank and residual degrees of freedom.", call. = FALSE)
  }
  design
}

#' @noRd
.validate_threads <- function(threads) {
  if (!is.numeric(threads) || length(threads) != 1L || !is.finite(threads) ||
        threads < 1 || threads != floor(threads)) {
    stop("'threads' must be a positive integer.", call. = FALSE)
  }
}

#' @noRd
.run_tool <- function(tool, args, stdout = TRUE) {
  check_tool(tool)
  diagnostics <- tempfile("bambamR-stderr-")
  on.exit(unlink(diagnostics), add = TRUE)
  result <- suppressWarnings(system2(tool, args = shQuote(as.character(args)),
                                     stdout = stdout, stderr = diagnostics))
  status <- if (is.character(result)) attr(result, "status") else result
  if (is.null(status)) status <- 0L
  if (length(status) != 1L || is.na(status) || status != 0L) {
    detail <- if (file.exists(diagnostics)) readLines(diagnostics, warn = FALSE) else character()
    stop(tool, " failed (exit ", status, "): ",
         paste(utils::tail(detail, 20L), collapse = "\n"), call. = FALSE)
  }
  result
}

#' @noRd
.require_product <- function(path) {
  if (!file.exists(path) || is.na(file.info(path)$size) || file.info(path)$size == 0) {
    stop("External command did not produce a non-empty output: ", path, call. = FALSE)
  }
}


# Internal helper: validate DE result
# ------------------------------------------------------------------
#' @noRd
validate_de_result <- function(de_result) {
  required <- c("gene", "log2fc", "pvalue", "padj")
  missing <- setdiff(required, colnames(de_result))
  if (length(missing) > 0L) {
    stop(
      "DE result missing required columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(de_result)
}


# Internal helper: standardize DE result column names
# ------------------------------------------------------------------
#' @noRd
standardize_de_result <- function(df, gene_col = "gene", log2fc_col = "log2fc",
                                  pvalue_col = "pvalue", padj_col = "padj",
                                  basemean_col = NULL) {
  result <- data.frame(
    gene    = df[[gene_col]],
    log2fc  = as.numeric(df[[log2fc_col]]),
    pvalue  = as.numeric(df[[pvalue_col]]),
    padj    = as.numeric(df[[padj_col]]),
    stringsAsFactors = FALSE
  )
  if (!is.null(basemean_col) && basemean_col %in% colnames(df)) {
    result$basemean <- as.numeric(df[[basemean_col]])
  }
  result
}


# S3 class: bb_result
# ------------------------------------------------------------------
#' Constructor for bb_result
#' @noRd
new_bb_result <- function(counts = NULL, metadata = NULL,
                          de_results = NULL, plots = list(),
                          qc = NULL, alignment_stats = NULL) {
  structure(
    list(
      counts          = counts,
      metadata        = metadata,
      de_results      = de_results,
      plots           = plots,
      qc              = qc,
      alignment_stats = alignment_stats
    ),
    class = "bb_result"
  )
}

#' Print method for bb_result
#' @param x A bb_result object.
#' @param ... Additional arguments (ignored).
#' @return `x`, invisibly; called for its side effect of printing a summary.
#' @export
print.bb_result <- function(x, ...) {
  cat("bambamR Result\n")
  cat("==============\n")
  if (!is.null(x$counts)) {
    cat("Genes:   ", nrow(x$counts), "\n")
    cat("Samples: ", ncol(x$counts), "\n")
  }
  if (!is.null(x$de_results)) {
    cat("DE results: available (",
        sum(x$de_results$padj < 0.05, na.rm = TRUE),
        " significant at FDR < 0.05)\n")
  }
  if (length(x$plots) > 0L) {
    cat("Plots: ", paste(names(x$plots), collapse = ", "), "\n")
  }
  invisible(x)
}
