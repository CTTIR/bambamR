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
#>    gene      log2fc    pvalue     padj  basemean
#> 1 gene1  0.12495253 0.4573411 0.992101  96.46432
#> 2 gene2  0.12823113 0.4336525 0.992101  97.29936
#> 3 gene3 -0.06742123 0.6617299 0.992101 100.65546
#> 4 gene4  0.06967958 0.6798876 0.992101  95.66183
#> 5 gene5 -0.07691483 0.6156012 0.992101 100.67061
#> 6 gene6 -0.11638884 0.4447977 0.992101  99.33898
# }
```
