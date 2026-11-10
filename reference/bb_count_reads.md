# Count Reads per Gene/Feature

Counts reads overlapping genomic features. Uses
[`GenomicAlignments::summarizeOverlaps()`](https://rdrr.io/pkg/GenomicAlignments/man/summarizeOverlaps-methods.html)
if available, otherwise shells out to `featureCounts` from the Subread
package.

## Usage

``` r
bb_count_reads(
  bam_paths,
  annotation,
  method = c("auto", "internal", "featureCounts"),
  threads = 4L,
  feature_type = "exon",
  attr_type = "gene_id",
  paired = FALSE
)
```

## Arguments

- bam_paths:

  Character vector. Paths to BAM files. Optional vector names supply
  unique sample IDs when BAM basenames are not unique.

- annotation:

  Character. Path to GTF/GFF annotation file.

- method:

  Character. `"auto"` tries GenomicAlignments first, then featureCounts;
  `"internal"` forces GenomicAlignments; `"featureCounts"` forces the
  external tool.

- threads:

  Integer. Number of threads for featureCounts. Default `4`.

- feature_type:

  Character. Feature type to count (GTF column 3). Default `"exon"`.

- attr_type:

  Character. Attribute to group by. Default `"gene_id"`.

- paired:

  Logical. Count paired-end fragments. Default `FALSE`. The
  featureCounts backend requires Subread 2.0.2 or later in paired mode.

## Value

A numeric matrix of counts (genes x samples). Rownames are gene IDs,
colnames are derived from BAM file names.

## Examples

``` r
# \donttest{
# Internal counting needs the optional Bioconductor packages Rsamtools,
# GenomicAlignments, GenomicFeatures and txdbmaker.
needed <- c("Rsamtools", "GenomicAlignments", "GenomicRanges",
            "GenomicFeatures", "txdbmaker", "SummarizedExperiment")
if (all(vapply(needed, requireNamespace, logical(1), quietly = TRUE))) {
  dir <- tempfile("bambamR-count-")
  dir.create(dir)
  sam <- file.path(dir, "s1.sam")
  writeLines(c(
    "@HD\tVN:1.6\tSO:coordinate", "@SQ\tSN:chr1\tLN:2000",
    "r1\t0\tchr1\t100\t60\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII",
    "r2\t0\tchr1\t300\t60\t8M\t*\t0\t0\tTTGCAACG\tIIIIIIII"
  ), sam)
  bam <- Rsamtools::asBam(sam, destination = file.path(dir, "s1"))
  gtf <- file.path(dir, "genes.gtf")
  writeLines(c(
    paste("chr1", "example", "exon", "90", "150", ".", "+", ".",
          "gene_id \"G1\"; transcript_id \"T1\";", sep = "\t"),
    paste("chr1", "example", "exon", "280", "350", ".", "+", ".",
          "gene_id \"G2\"; transcript_id \"T2\";", sep = "\t")
  ), gtf)
  bb_count_reads(bam, gtf, method = "internal")
  unlink(dir, recursive = TRUE)
}
#> Counting reads with GenomicAlignments...
#> Import genomic features from the file as a GRanges object ... 
#> OK
#> Prepare the 'metadata' data frame ... 
#> OK
#> Make the TxDb object ... 
#> Warning: genome version information is not available for this TxDb object
#> OK
# }
```
