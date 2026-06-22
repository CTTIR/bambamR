# Export DE Results to CSV

Writes a differential expression result data.frame to a CSV file.

## Usage

``` r
bb_export_csv(de_result, path)
```

## Arguments

- de_result:

  A data.frame with DE results (columns: gene, log2fc, pvalue, padj, and
  optionally basemean).

- path:

  Character. File path for the output `.csv` file.

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
tmp <- tempfile(fileext = ".csv")
bb_export_csv(de, tmp)
#> Saved to: /tmp/RtmpfRpUwx/file1a84133150cf.csv
```
