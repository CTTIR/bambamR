# Generate QC Metrics

Computes quality control metrics from FASTQ and/or BAM files.

## Usage

``` r
bb_qc(fastq_path = NULL, bam_path = NULL, use_fastqc = TRUE)
```

## Arguments

- fastq_path:

  Character or NULL. Path(s) to FASTQ file(s).

- bam_path:

  Character or NULL. Path(s) to BAM file(s).

- use_fastqc:

  Logical. Try to use system FastQC if available. Default `TRUE`.

## Value

A `bb_qc` object containing:

- read_counts:

  Total reads per file

- quality_scores:

  Per-position quality score summary

- gc_content:

  GC content distribution

- read_lengths:

  Read length distribution

- mapping_rate:

  Mapping rate from BAM files (if provided)

## Examples

``` r
# \donttest{
# Create a test FASTQ
tmp <- tempfile(fileext = ".fastq")
writeLines(c(
  "@read1", "ACGTACGT", "+", "IIIIIIII",
  "@read2", "GCGCGCGC", "+", "HHHHHHHH"
), tmp)
qc <- bb_qc(fastq_path = tmp)
qc
#> bambamR QC Summary
#> ==================
#> Files analyzed: 1 
#>   file334b178a0d46.fastq: 2 reads
# }
```
