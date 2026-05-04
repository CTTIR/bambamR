#' Differential Expression with DESeq2
#'
#' Runs DESeq2 differential expression analysis. Requires the `DESeq2` package.
#'
#' @param counts Numeric matrix. Raw (unnormalized) count matrix.
#' @param coldata A data.frame with sample information. Rownames must match
#'   column names of `counts`. Must contain the variables referenced in
#'   `design`.
#' @param design A formula specifying the design. Default `~ condition`.
#' @param contrast Character vector of length 3: `c("variable", "numerator",
#'   "denominator")`. If `NULL`, the last variable in the design is used.
#' @param alpha Numeric. FDR threshold for independent filtering.
#'   Default `0.05`.
#'
#' @return A data.frame with columns `gene`, `log2fc`, `pvalue`, `padj`,
#'   `basemean`.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("DESeq2", quietly = TRUE)) {
#'   counts <- matrix(rpois(600, 100), nrow = 100, ncol = 6,
#'     dimnames = list(paste0("gene", 1:100), paste0("S", 1:6)))
#'   coldata <- data.frame(
#'     condition = factor(rep(c("ctrl", "treat"), each = 3)),
#'     row.names = paste0("S", 1:6)
#'   )
#'   result <- bb_deseq2(counts, coldata)
#'   head(result)
#' }
#' }
#'
#' @export
bb_deseq2 <- function(counts, coldata, design = ~ condition,
                       contrast = NULL, alpha = 0.05) {

  check_pkg("DESeq2", reason = "for DESeq2 differential expression")
  counts <- validate_counts(counts)

  # Ensure matching samples
  shared <- intersect(colnames(counts), rownames(coldata))
  if (length(shared) < 2L) {
    stop("Need at least 2 matching samples between counts and coldata.",
         call. = FALSE)
  }
  counts <- counts[, shared, drop = FALSE]
  coldata <- coldata[shared, , drop = FALSE]

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = counts,
    colData   = coldata,
    design    = design
  )
  dds <- DESeq2::DESeq(dds)

  if (!is.null(contrast)) {
    res <- DESeq2::results(dds, contrast = contrast, alpha = alpha)
  } else {
    res <- DESeq2::results(dds, alpha = alpha)
  }

  res_df <- as.data.frame(res)
  standardize_de_result(
    data.frame(
      gene     = rownames(res_df),
      log2fc   = res_df$log2FoldChange,
      pvalue   = res_df$pvalue,
      padj     = res_df$padj,
      basemean = res_df$baseMean,
      stringsAsFactors = FALSE
    ),
    gene_col     = "gene",
    log2fc_col   = "log2fc",
    pvalue_col   = "pvalue",
    padj_col     = "padj",
    basemean_col = "basemean"
  )
}


#' Differential Expression with edgeR
#'
#' Runs edgeR quasi-likelihood differential expression analysis. Requires
#' the `edgeR` package.
#'
#' @param counts Numeric matrix. Raw count matrix.
#' @param group Factor or character vector. Group assignment for each sample.
#' @param design A design matrix. If `NULL`, a simple group design is created.
#'
#' @return A data.frame with columns `gene`, `log2fc`, `pvalue`, `padj`.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("edgeR", quietly = TRUE)) {
#'   counts <- matrix(rpois(600, 100), nrow = 100, ncol = 6,
#'     dimnames = list(paste0("gene", 1:100), paste0("S", 1:6)))
#'   group <- factor(rep(c("ctrl", "treat"), each = 3))
#'   result <- bb_edger(counts, group)
#'   head(result)
#' }
#' }
#'
#' @export
bb_edger <- function(counts, group, design = NULL) {

  check_pkg("edgeR", reason = "for edgeR differential expression")
  counts <- validate_counts(counts)

  if (length(group) != ncol(counts)) {
    stop("Length of 'group' must equal number of columns in 'counts'.",
         call. = FALSE)
  }

  group <- factor(group)
  dge <- edgeR::DGEList(counts = counts, group = group)
  dge <- edgeR::calcNormFactors(dge)

  if (is.null(design)) {
    design <- stats::model.matrix(~ group)
  }

  dge <- edgeR::estimateDisp(dge, design)
  fit <- edgeR::glmQLFit(dge, design)
  qlf <- edgeR::glmQLFTest(fit, coef = ncol(design))
  tt <- edgeR::topTags(qlf, n = nrow(counts), sort.by = "none")$table

  standardize_de_result(
    data.frame(
      gene   = rownames(tt),
      log2fc = tt$logFC,
      pvalue = tt$PValue,
      padj   = tt$FDR,
      stringsAsFactors = FALSE
    )
  )
}


#' Differential Expression with limma-voom
#'
#' Runs limma-voom differential expression analysis. Requires the `limma`
#' and `edgeR` packages.
#'
#' @param counts Numeric matrix. Raw count matrix.
#' @param design A design matrix (e.g., from [stats::model.matrix()]).
#' @param contrast A contrast vector or matrix. If `NULL`, the last
#'   coefficient is tested.
#'
#' @return A data.frame with columns `gene`, `log2fc`, `pvalue`, `padj`.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("limma", quietly = TRUE) &&
#'     requireNamespace("edgeR", quietly = TRUE)) {
#'   counts <- matrix(rpois(600, 100), nrow = 100, ncol = 6,
#'     dimnames = list(paste0("gene", 1:100), paste0("S", 1:6)))
#'   group <- factor(rep(c("ctrl", "treat"), each = 3))
#'   design <- model.matrix(~ group)
#'   result <- bb_limma_voom(counts, design)
#'   head(result)
#' }
#' }
#'
#' @export
bb_limma_voom <- function(counts, design, contrast = NULL) {

  check_pkg("limma", reason = "for limma-voom analysis")
  check_pkg("edgeR", reason = "for voom normalization")
  counts <- validate_counts(counts)

  dge <- edgeR::DGEList(counts = counts)
  dge <- edgeR::calcNormFactors(dge)

  v <- limma::voom(dge, design, plot = FALSE)
  fit <- limma::lmFit(v, design)

  if (!is.null(contrast)) {
    fit <- limma::contrasts.fit(fit, contrasts = contrast)
  }

  fit <- limma::eBayes(fit)
  tt <- limma::topTable(fit, coef = ncol(design), number = nrow(counts),
                        sort.by = "none")

  standardize_de_result(
    data.frame(
      gene   = rownames(tt),
      log2fc = tt$logFC,
      pvalue = tt$P.Value,
      padj   = tt$adj.P.Val,
      stringsAsFactors = FALSE
    )
  )
}
