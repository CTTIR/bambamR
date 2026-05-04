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
                 "ComplexHeatmap")
  available <- vapply(bioc_pkgs, requireNamespace, logical(1), quietly = TRUE)

  if (all(available)) {
    packageStartupMessage("bambamR: Full mode (all Bioconductor packages available)")
  } else {
    missing <- bioc_pkgs[!available]
    packageStartupMessage(
      "bambamR: Minimal mode (CRAN-only). ",
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
                   "Biobase", "ComplexHeatmap")
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
  counts
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
