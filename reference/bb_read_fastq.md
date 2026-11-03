# Read FASTQ Files

Reads FASTQ files using ShortRead (if available) or a minimal base-R
parser as fallback. Supports `.gz` compressed files.

## Usage

``` r
bb_read_fastq(path, n = NULL)
```

## Arguments

- path:

  Character. Path to a FASTQ or FASTQ.gz file.

- n:

  Integer or NULL. Number of initial reads to read. `NULL` reads all. A
  finite limit uses the four-line base parser without reading the full
  file.

## Value

A data.frame with columns:

- id:

  Read identifier

- sequence:

  Nucleotide sequence

- quality:

  Quality string (Phred+33 encoded)

## Examples

``` r
# \donttest{
# Create a temporary FASTQ file
tmp <- tempfile(fileext = ".fastq")
writeLines(c(
  "@read1", "ACGTACGT", "+", "IIIIIIII",
  "@read2", "TGCATGCA", "+", "HHHHHHHH"
), tmp)
reads <- bb_read_fastq(tmp)
reads
#>      id sequence  quality
#> 1 read1 ACGTACGT IIIIIIII
#> 2 read2 TGCATGCA HHHHHHHH
# }
```
