# Export DE Results to TSV

Writes a differential expression result data.frame to a TSV file.

## Usage

``` r
bb_export_tsv(de_result, path)
```

## Arguments

- de_result:

  A data.frame with DE results.

- path:

  Character. File path for the output `.tsv` file.

## Value

The `path` invisibly.

## Examples

``` r
de <- data.frame(
  gene = paste0("gene", 1:5),
  log2fc = rnorm(5),
  pvalue = runif(5, 0, 0.1),
  padj = runif(5, 0, 0.1)
)
tmp <- tempfile(fileext = ".tsv")
bb_export_tsv(de, tmp)
#> Saved to: /tmp/RtmpzM8Az8/file334ba8956fe.tsv
```
