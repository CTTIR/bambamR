# Differential Expression with edgeR

Runs edgeR quasi-likelihood differential expression analysis. Requires
the `edgeR` package.

## Usage

``` r
bb_edger(counts, group, design = NULL)
```

## Arguments

- counts:

  Numeric matrix. Raw count matrix.

- group:

  Factor or character vector. Group assignment for each sample.

- design:

  A design matrix. If `NULL`, a simple group design is created.

## Value

A data.frame with columns `gene`, `log2fc`, `pvalue`, `padj`.

## Examples

``` r
# \donttest{
if (requireNamespace("edgeR", quietly = TRUE)) {
  counts <- matrix(rpois(600, 100), nrow = 100, ncol = 6,
    dimnames = list(paste0("gene", 1:100), paste0("S", 1:6)))
  group <- factor(rep(c("ctrl", "treat"), each = 3))
  result <- bb_edger(counts, group)
  head(result)
}
#> calcNormFactors has been renamed to normLibSizes
#>    gene      log2fc    pvalue      padj
#> 1 gene1  0.03672190 0.7383797 0.9569488
#> 2 gene2 -0.03460213 0.7565203 0.9569488
#> 3 gene3  0.05702904 0.6369530 0.9368944
#> 4 gene4 -0.05426078 0.6555098 0.9368944
#> 5 gene5  0.16545800 0.1808186 0.8210933
#> 6 gene6 -0.06047316 0.5845282 0.9368944
# }
```
