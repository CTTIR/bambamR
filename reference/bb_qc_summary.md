# QC Summary

Returns a summary data.frame of QC metrics.

## Usage

``` r
bb_qc_summary(qc_object)
```

## Arguments

- qc_object:

  A `bb_qc` object from
  [`bb_qc()`](https://cttir.github.io/bambamR/reference/bb_qc.md).

## Value

A data.frame with one row per file.

## Examples

``` r
# \donttest{
tmp <- tempfile(fileext = ".fastq")
writeLines(c("@r1", "ACGT", "+", "IIII"), tmp)
qc <- bb_qc(fastq_path = tmp)
bb_qc_summary(qc)
#>                     file total_reads median_gc mapping_rate
#> 1 file33456c2f42ca.fastq           1       0.5           NA
# }
```
