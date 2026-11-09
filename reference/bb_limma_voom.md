# Differential Expression with limma-voom

Runs limma-voom differential expression analysis. Requires the `limma`
and `edgeR` packages.

## Usage

``` r
bb_limma_voom(counts, design, contrast = NULL)
```

## Arguments

- counts:

  Numeric matrix. Raw count matrix.

- design:

  A design matrix (e.g., from
  [`stats::model.matrix()`](https://rdrr.io/r/stats/model.matrix.html)).

- contrast:

  A contrast vector or one-column matrix. If `NULL`, the last design
  coefficient is tested. One contrast is returned per call.

## Value

A data.frame with columns `gene`, `log2fc`, `pvalue`, `padj`.

## Examples

``` r
# \donttest{
if (requireNamespace("limma", quietly = TRUE) &&
    requireNamespace("edgeR", quietly = TRUE)) {
  counts <- matrix(rpois(600, 100), nrow = 100, ncol = 6,
    dimnames = list(paste0("gene", 1:100), paste0("S", 1:6)))
  group <- factor(rep(c("ctrl", "treat"), each = 3))
  design <- model.matrix(~ group)
  result <- bb_limma_voom(counts, design)
  head(result)
}
#> calcNormFactors has been renamed to normLibSizes
#>    gene      log2fc     pvalue      padj
#> 1 gene1 -0.07229408 0.52407853 0.9092510
#> 2 gene2 -0.04190298 0.71594977 0.9549223
#> 3 gene3  0.04908073 0.60525976 0.9315685
#> 4 gene4  0.08188223 0.56818642 0.9164297
#> 5 gene5 -0.05009243 0.60471181 0.9315685
#> 6 gene6  0.22960595 0.02025713 0.6498947
# }
```
