# Run Full bambamR Pipeline

Orchestrates the complete RNA-seq analysis pipeline from FASTQ/BAM files
or a count matrix through alignment, counting, normalization,
differential expression, and visualization.

## Usage

``` r
bb_pipeline(
  fastq_dir = NULL,
  bam_dir = NULL,
  count_matrix = NULL,
  output_dir = tempfile("bambamR_output"),
  genome_index = NULL,
  annotation = NULL,
  sample_info = NULL,
  aligner = "STAR",
  de_method = c("DESeq2", "edgeR", "limma"),
  design = ~condition,
  skip = character(0),
  threads = 4L
)
```

## Arguments

- fastq_dir:

  Character or NULL. Directory containing FASTQ files. If provided, the
  pipeline starts from alignment. Each file is one single-end sample,
  named by removing its FASTQ extension. For paired-end reads use
  [`bb_align()`](https://cttir.github.io/bambamR/reference/bb_align.md)
  with explicit pairs, then
  [`bb_count_reads()`](https://cttir.github.io/bambamR/reference/bb_count_reads.md).

- bam_dir:

  Character or NULL. Directory containing BAM files. If provided, the
  pipeline starts from read counting.

- count_matrix:

  A numeric matrix or NULL. Pre-computed count matrix. If provided, the
  pipeline starts from normalization/DE.

- output_dir:

  Character. Output directory for results. Defaults to a new directory
  below [`tempdir()`](https://rdrr.io/r/base/tempfile.html); supply a
  permanent path to keep the results. Existing result files are never
  overwritten.

- genome_index:

  Character or NULL. Path to genome index. Required if starting from
  FASTQ.

- annotation:

  Character or NULL. Path to GTF annotation. Required if starting from
  FASTQ or BAM.

- sample_info:

  A data.frame with sample metadata. Must contain a `condition` column
  (or the variable in `design`). Rownames should match sample
  identifiers.

- aligner:

  Character. Aligner for FASTQ alignment. Default `"STAR"`.

- de_method:

  Character. DE method: `"DESeq2"`, `"edgeR"`, or `"limma"`. Default
  `"DESeq2"`.

- design:

  A formula for DE analysis. Default `~ condition`.

- skip:

  Character vector. Steps to skip: `"qc"`, `"align"`, `"count"`, `"de"`,
  `"viz"`. Default `character(0)`.

- threads:

  Integer. Number of threads. Default `4`.

## Value

A `bb_result` object containing counts, metadata, DE results, and plots.

## Details

The pipeline allows entry at any stage:

- **From FASTQ**: requires `fastq_dir`, `genome_index`, `annotation`,
  and `sample_info`

- **From BAM**: requires `bam_dir`, `annotation`, and `sample_info`

- **From counts**: requires `count_matrix` and `sample_info`

## Examples

``` r
# Starting from a count matrix; skipping differential expression keeps
# this example free of optional Bioconductor packages.
set.seed(42)
counts <- matrix(rpois(600, 100), nrow = 100, ncol = 6,
  dimnames = list(paste0("gene", 1:100), paste0("S", 1:6)))
sample_info <- data.frame(
  condition = factor(rep(c("ctrl", "treat"), each = 3)),
  row.names = paste0("S", 1:6)
)
out <- file.path(tempdir(), "bambamR-pipeline-example")
result <- bb_pipeline(count_matrix = counts, sample_info = sample_info,
                      output_dir = out, skip = "de")
#> === Step 4: Normalization ===
#>   CPM normalization complete
#> === Step 6: Generating Plots ===
#>   Generated: pca, heatmap
#> Saved to: /tmp/RtmpPdJiwa/bambamR-pipeline-example/bambamR_result.rds
#> 
#> === Pipeline complete ===
result
#> bambamR Result
#> ==============
#> Genes:    100 
#> Samples:  6 
#> Plots:  pca, heatmap 
unlink(out, recursive = TRUE)

# \donttest{
# Full run including differential expression with edgeR.
if (requireNamespace("edgeR", quietly = TRUE)) {
  out <- file.path(tempdir(), "bambamR-pipeline-edger")
  result <- bb_pipeline(count_matrix = counts, sample_info = sample_info,
                        output_dir = out, de_method = "edgeR")
  head(result$de_results)
  unlink(out, recursive = TRUE)
}
#> === Step 4: Normalization ===
#>   CPM normalization complete
#> === Step 5: Differential Expression (edgeR) ===
#> calcNormFactors has been renamed to normLibSizes
#>   Found 0 significant genes (FDR < 0.05)
#> === Step 6: Generating Plots ===
#>   Generated: pca, heatmap, volcano
#> Saved to: /tmp/RtmpPdJiwa/bambamR-pipeline-edger/bambamR_result.rds
#> 
#> === Pipeline complete ===
# }
```
