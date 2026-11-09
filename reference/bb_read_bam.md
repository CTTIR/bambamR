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
# Requires a BAM file
# bam_df <- bb_read_bam("aligned.bam")
# }
```
