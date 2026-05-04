#' Normalize Count Matrix
#'
#' Supports CPM, TPM, TMM (edgeR), and RLE (DESeq2) normalization.
#' CPM and TPM are always available. TMM requires edgeR, RLE requires DESeq2.
#'
#' @param counts Numeric matrix. Raw count matrix with gene IDs as rownames
#'   and sample names as colnames.
#' @param method Character. One of `"cpm"`, `"tpm"`, `"tmm"`, `"rle"`.
#' @param gene_lengths Numeric vector. Gene lengths in base pairs. Required
#'   for TPM normalization. Must be named or in the same order as rownames
#'   of `counts`.
#'
#' @return A normalized numeric matrix with the same dimensions as `counts`.
#'
#' @details
#' - **CPM** (Counts Per Million): `counts / library_size * 1e6`. Pure base R.
#' - **TPM** (Transcripts Per Million): Normalizes by gene length then library
#'   size. Pure base R but requires `gene_lengths`.
#' - **TMM** (Trimmed Mean of M-values): Uses `edgeR::calcNormFactors()`.
#'   Requires the `edgeR` package.
#' - **RLE** (Relative Log Expression): Uses `DESeq2::estimateSizeFactors()`.
#'   Requires the `DESeq2` package.
#'
#' @examples
#' counts <- matrix(
#'   rpois(600, lambda = 100),
#'   nrow = 100, ncol = 6,
#'   dimnames = list(paste0("gene", 1:100), paste0("sample", 1:6))
#' )
#' cpm <- bb_normalize(counts, method = "cpm")
#' head(cpm)
#'
#' # TPM requires gene lengths
#' gene_lengths <- sample(500:5000, 100)
#' tpm <- bb_normalize(counts, method = "tpm", gene_lengths = gene_lengths)
#'
#' @export
bb_normalize <- function(counts, method = c("cpm", "tpm", "tmm", "rle"),
                         gene_lengths = NULL) {

  method <- match.arg(method)
  counts <- validate_counts(counts)

  switch(method,
    cpm = .normalize_cpm(counts),
    tpm = .normalize_tpm(counts, gene_lengths),
    tmm = .normalize_tmm(counts),
    rle = .normalize_rle(counts)
  )
}


# CPM: pure base R
# ------------------------------------------------------------------
#' @noRd
.normalize_cpm <- function(counts) {
  lib_sizes <- colSums(counts)
  t(t(counts) / lib_sizes) * 1e6
}


# TPM: pure base R, needs gene_lengths
# ------------------------------------------------------------------
#' @noRd
.normalize_tpm <- function(counts, gene_lengths) {
  if (is.null(gene_lengths)) {
    stop("'gene_lengths' is required for TPM normalization.", call. = FALSE)
  }
  if (length(gene_lengths) != nrow(counts)) {
    stop("Length of 'gene_lengths' must match number of rows in 'counts'.",
         call. = FALSE)
  }
  # Reads per kilobase
  rpk <- counts / (gene_lengths / 1000)
  # Per-million scaling factor
  scaling <- colSums(rpk) / 1e6
  t(t(rpk) / scaling)
}


# TMM: requires edgeR
# ------------------------------------------------------------------
#' @noRd
.normalize_tmm <- function(counts) {
  check_pkg("edgeR", reason = "for TMM normalization")
  dge <- edgeR::DGEList(counts = counts)

  dge <- edgeR::calcNormFactors(dge, method = "TMM")
  edgeR::cpm(dge)
}


# RLE: requires DESeq2
# ------------------------------------------------------------------
#' @noRd
.normalize_rle <- function(counts) {
  check_pkg("DESeq2", reason = "for RLE normalization")
  # Create a minimal DESeqDataSet
  coldata <- data.frame(condition = rep("A", ncol(counts)),
                        row.names = colnames(counts))
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = counts,
    colData = coldata,
    design = ~ 1
  )
  dds <- DESeq2::estimateSizeFactors(dds)
  DESeq2::counts(dds, normalized = TRUE)
}
