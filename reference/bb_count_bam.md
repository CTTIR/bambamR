# Count Reads in BAM File

Returns the total number of reads in a BAM file using Rsamtools or
system samtools.

## Usage

``` r
bb_count_bam(path)
```

## Arguments

- path:

  Character. Path to a BAM file.

## Value

A numeric value: total number of reads.

## Examples

``` r
# \donttest{
# Building a BAM file requires the optional Bioconductor package Rsamtools,
# whose loading can take several seconds.
if (requireNamespace("Rsamtools", quietly = TRUE)) {
  dir <- tempfile("bambamR-bam-")
  dir.create(dir)
  sam <- file.path(dir, "reads.sam")
  writeLines(c(
    "@HD\tVN:1.6\tSO:coordinate", "@SQ\tSN:chr1\tLN:2000",
    "r1\t0\tchr1\t100\t60\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII",
    "r2\t0\tchr1\t300\t60\t8M\t*\t0\t0\tTTGCAACG\tIIIIIIII"
  ), sam)
  bam <- Rsamtools::asBam(sam, destination = file.path(dir, "reads"))
  bb_count_bam(bam)
  unlink(dir, recursive = TRUE)
}
# }
```
