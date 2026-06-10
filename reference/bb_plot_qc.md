# QC Visualization Panel

Creates a panel of QC plots from a `bb_qc` object, including read
quality per position, GC content distribution, and read length
distribution.

## Usage

``` r
bb_plot_qc(qc_object, which = c("all", "quality", "gc", "length"))
```

## Arguments

- qc_object:

  A `bb_qc` object from
  [`bb_qc()`](https://cttir.github.io/bambamR/reference/bb_qc.md).

- which:

  Character. Which plot to generate: `"quality"`, `"gc"`, `"length"`, or
  `"all"`. Default `"all"`.

## Value

A [ggplot2::ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html)
object (or a patchwork composition if `which = "all"` and patchwork is
available).

## Examples

``` r
# \donttest{
tmp <- tempfile(fileext = ".fastq")
writeLines(c(
  "@r1", "ACGTACGT", "+", "IIIIIIII",
  "@r2", "GCGCGCGC", "+", "HHHHHHHH",
  "@r3", "ATATGCGC", "+", "GGGGFFFF"
), tmp)
qc <- bb_qc(fastq_path = tmp)
bb_plot_qc(qc)

# }
```
