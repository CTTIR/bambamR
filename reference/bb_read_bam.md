# Read BAM Files

Reads BAM files using Rsamtools (if available) or the system `samtools`
command as fallback.

## Usage

``` r
bb_read_bam(
  path,
  index = TRUE,
  what = c("qname", "flag", "rname", "pos", "mapq", "cigar")
)
```

## Arguments

- path:

  Character. Path to a BAM file.

- index:

  Logical. Whether to use/create a BAM index. Default `TRUE`.

- what:

  Character vector. Fields to extract. Default includes basic alignment
  fields.

## Value

A data.frame with alignment information including columns: `qname`,
`flag`, `rname`, `pos`, `mapq`, `cigar`.

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
  bb_read_bam(bam)
  unlink(dir, recursive = TRUE)
}
# }
```
