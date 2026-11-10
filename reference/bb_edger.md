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
#>    gene      log2fc     pvalue      padj
#> 1 gene1  0.14622487 0.21590888 0.8996204
#> 2 gene2 -0.15832656 0.17448249 0.8996204
#> 3 gene3 -0.03622465 0.76026904 0.9711239
#> 4 gene4  0.21126922 0.08146277 0.8291105
#> 5 gene5  0.08666288 0.51657260 0.9415020
#> 6 gene6  0.10569611 0.34988150 0.9259242
# }
```
