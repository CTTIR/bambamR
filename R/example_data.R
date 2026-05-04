#' Load Example Count Matrix
#'
#' Loads a synthetic RNA-seq count matrix (200 genes x 10 samples) with
#' sample metadata. The first 20 genes have simulated differential expression
#' between control and treatment groups.
#'
#' @return A list with components:
#' \describe{
#'   \item{counts}{A 200 x 10 integer matrix of raw counts. Rownames are gene
#'     symbols, colnames are sample IDs.}
#'   \item{metadata}{A data.frame with columns `condition`, `batch`, `sex`,
#'     `age`. Rownames match colnames of `counts`.}
#' }
#'
#' @examples
#' ex <- bb_example_counts()
#' dim(ex$counts)
#' head(ex$metadata)
#'
#' # Quick pipeline from example data
#' cpm <- bb_normalize(ex$counts, method = "cpm")
#' bb_pca(cpm, ex$metadata, color_by = "condition")
#'
#' @export
bb_example_counts <- function() {
  path <- system.file("extdata", "example_counts.rds", package = "bambamR")
  if (path == "") {
    stop("Example data not found. Try reinstalling bambamR.", call. = FALSE)
  }
  readRDS(path)
}


#' Load Example Mutation Data
#'
#' Loads synthetic mutation data (300 mutations across 50 samples and 20
#' cancer-related genes) suitable for oncoplot visualization, along with
#' clinical annotation data.
#'
#' @return A list with components:
#' \describe{
#'   \item{mutations}{A data.frame with columns `sample`, `gene`,
#'     `mutation_type`.}
#'   \item{clinical}{A data.frame with columns `Stage`, `Gender`, `Smoking`.
#'     Rownames are sample IDs.}
#' }
#'
#' @examples
#' ex <- bb_example_mutations()
#' head(ex$mutations)
#'
#' # Basic oncoplot
#' bb_oncoplot(ex$mutations, n_genes = 10)
#'
#' # With clinical annotations
#' bb_oncoplot(ex$mutations, n_genes = 10, annotation_df = ex$clinical)
#'
#' @export
bb_example_mutations <- function() {
  path <- system.file("extdata", "example_mutations.rds", package = "bambamR")
  if (path == "") {
    stop("Example data not found. Try reinstalling bambamR.", call. = FALSE)
  }
  readRDS(path)
}


#' Load Example DE Results
#'
#' Loads pre-computed differential expression results (500 genes) with
#' realistic log2 fold-changes and p-values. The first 20 genes are
#' simulated as significantly differentially expressed.
#'
#' @return A data.frame with columns `gene`, `log2fc`, `pvalue`, `padj`,
#'   `basemean`.
#'
#' @examples
#' de <- bb_example_de()
#' head(de[order(de$padj), ])
#'
#' # Volcano plot
#' bb_volcano(de)
#'
#' # MA plot
#' bb_ma_plot(de)
#'
#' @export
bb_example_de <- function() {
  path <- system.file("extdata", "example_de_results.rds", package = "bambamR")
  if (path == "") {
    stop("Example data not found. Try reinstalling bambamR.", call. = FALSE)
  }
  readRDS(path)
}
