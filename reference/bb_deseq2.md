# Differential Expression with DESeq2

Runs DESeq2 differential expression analysis. Requires the `DESeq2`
package.

## Usage

``` r
bb_deseq2(counts, coldata, design = ~condition, contrast = NULL, alpha = 0.05)
```

## Arguments

- counts:

  Numeric matrix. Raw (unnormalized) count matrix.

- coldata:

  A data.frame with sample information. Rownames must match column names
  of `counts`. Must contain the variables referenced in `design`.

- design:

  A formula specifying the design. Default `~ condition`.

- contrast:

  Character vector of length 3:
  `c("variable", "numerator", "denominator")`. If `NULL`, the last
  variable in the design is used.

- alpha:

  Numeric. FDR threshold for independent filtering. Default `0.05`.

## Value

A data.frame with columns `gene`, `log2fc`, `pvalue`, `padj`,
`basemean`.

## Examples

``` r
# \donttest{
if (requireNamespace("DESeq2", quietly = TRUE)) {
  counts <- matrix(rpois(600, 100), nrow = 100, ncol = 6,
    dimnames = list(paste0("gene", 1:100), paste0("S", 1:6)))
  coldata <- data.frame(
    condition = factor(rep(c("ctrl", "treat"), each = 3)),
    row.names = paste0("S", 1:6)
  )
  result <- bb_deseq2(counts, coldata)
  head(result)
}
#> estimating size factors
#> estimating dispersions
#> gene-wise dispersion estimates
#> mean-dispersion relationship
#> -- note: fitType='parametric', but the dispersion trend was not well captured by the
#>    function: y = a/x + b, and a local regression fit was automatically substituted.
#>    specify fitType='local' or 'mean' to avoid this message next time.
#> final dispersion estimates
#> fitting model and testing
#>    gene     log2fc    pvalue      padj  basemean
#> 1 gene1 -0.1064135 0.5614828 0.9911937  94.46272
#> 2 gene2  0.1211631 0.4707732 0.9911937 102.23745
#> 3 gene3 -0.1815007 0.2568900 0.9911937  97.79708
#> 4 gene4 -0.1221147 0.4396177 0.9911937 102.58136
#> 5 gene5 -0.2027257 0.1972465 0.9911937 101.79267
#> 6 gene6 -0.1010836 0.5225444 0.9911937  89.81135
# }
```
