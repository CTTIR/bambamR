# Summary method for bb_qc

Summary method for bb_qc

## Usage

``` r
# S3 method for class 'bb_qc'
summary(object, ...)
```

## Arguments

- object:

  A bb_qc object.

- ...:

  Additional arguments (ignored).

## Value

A data frame with one row per analyzed file, as returned by
[`bb_qc_summary()`](https://cttir.github.io/bambamR/reference/bb_qc_summary.md).

## Examples

``` r
# \donttest{
# bb_qc() reads FASTQ files with the optional 'ShortRead' package when it
# is installed; loading it can take several seconds.
tmp <- tempfile(fileext = ".fastq")
writeLines(c("@r1", "ACGT", "+", "IIII"), tmp)
qc <- bb_qc(fastq_path = tmp, use_fastqc = FALSE)
summary(qc)
#>                     file total_reads median_gc mapping_rate
#> 1 file18956dbcc050.fastq           1       0.5           NA
unlink(tmp)
# }
```
