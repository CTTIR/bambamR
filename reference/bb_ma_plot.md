# MA Plot

Creates an MA plot (log2 fold-change vs. mean expression) from
differential expression results.

## Usage

``` r
bb_ma_plot(
  de_result,
  p_cutoff = 0.05,
  point_size = 1,
  colors = c(sig = "#D73027", ns = "grey70")
)
```

## Arguments

- de_result:

  A data.frame with columns `gene`, `log2fc`, `pvalue`, `padj`, and
  `basemean`.

- p_cutoff:

  Numeric. Adjusted p-value cutoff for coloring significant genes.
  Default `0.05`.

- point_size:

  Numeric. Size of points. Default `1`.

- colors:

  Named character vector with `"sig"` and `"ns"` colors.

## Value

A [ggplot2::ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html)
object.

## Examples

``` r
de <- data.frame(
  gene = paste0("gene", 1:200),
  log2fc = rnorm(200, 0, 2),
  pvalue = 10^(-runif(200, 0, 5)),
  padj = 10^(-runif(200, 0, 4)),
  basemean = 10^runif(200, 1, 4)
)
bb_ma_plot(de)

```
