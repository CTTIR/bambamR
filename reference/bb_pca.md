# PCA Plot

Creates a PCA plot from a normalized count matrix with sample metadata.

## Usage

``` r
bb_pca(
  counts,
  metadata,
  color_by,
  shape_by = NULL,
  n_genes = 500L,
  label = FALSE,
  point_size = 3
)
```

## Arguments

- counts:

  Numeric matrix. Normalized count matrix (genes x samples).

- metadata:

  A data.frame with sample information. Rownames must match column names
  of `counts`.

- color_by:

  Character. Column name in `metadata` to use for coloring points.

- shape_by:

  Character or NULL. Column name in `metadata` for point shapes.

- n_genes:

  Integer. Number of top variable genes to use for PCA. Default `500`.

- label:

  Logical. Whether to label sample points. Default `FALSE`.

- point_size:

  Numeric. Size of points. Default `3`.

## Value

A [ggplot2::ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html)
object.

## Examples

``` r
set.seed(42)
counts <- matrix(rpois(600, 100), nrow = 100, ncol = 6,
  dimnames = list(paste0("gene", 1:100), paste0("S", 1:6)))
meta <- data.frame(
  condition = rep(c("Control", "Treatment"), each = 3),
  row.names = paste0("S", 1:6)
)
bb_pca(counts, meta, color_by = "condition")

```
