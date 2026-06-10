# Changelog

## bambamR 0.1.0

### New features

- Initial release of bambamR.
- **Data import**:
  [`bb_read_fastq()`](https://cttir.github.io/bambamR/reference/bb_read_fastq.md)
  and
  [`bb_read_bam()`](https://cttir.github.io/bambamR/reference/bb_read_bam.md)
  with ShortRead/Rsamtools or base-R/system-samtools fallback.
- **Quality control**:
  [`bb_qc()`](https://cttir.github.io/bambamR/reference/bb_qc.md),
  [`bb_qc_summary()`](https://cttir.github.io/bambamR/reference/bb_qc_summary.md),
  and
  [`bb_plot_qc()`](https://cttir.github.io/bambamR/reference/bb_plot_qc.md)
  for per-read quality, GC content, and read-length distributions.
- **Alignment**:
  [`bb_align()`](https://cttir.github.io/bambamR/reference/bb_align.md)
  wrapping STAR, HISAT2, and minimap2.
- **Read counting**:
  [`bb_count_reads()`](https://cttir.github.io/bambamR/reference/bb_count_reads.md)
  via GenomicAlignments or featureCounts.
- **Normalization**:
  [`bb_normalize()`](https://cttir.github.io/bambamR/reference/bb_normalize.md)
  supporting CPM, TPM (base R), TMM (edgeR), and RLE (DESeq2).
- **Differential expression**:
  [`bb_deseq2()`](https://cttir.github.io/bambamR/reference/bb_deseq2.md),
  [`bb_edger()`](https://cttir.github.io/bambamR/reference/bb_edger.md),
  and
  [`bb_limma_voom()`](https://cttir.github.io/bambamR/reference/bb_limma_voom.md)
  returning standardized result data.frames.
- **Visualization**:
  [`bb_oncoplot()`](https://cttir.github.io/bambamR/reference/bb_oncoplot.md),
  [`bb_volcano()`](https://cttir.github.io/bambamR/reference/bb_volcano.md),
  [`bb_heatmap()`](https://cttir.github.io/bambamR/reference/bb_heatmap.md),
  [`bb_pca()`](https://cttir.github.io/bambamR/reference/bb_pca.md),
  [`bb_ma_plot()`](https://cttir.github.io/bambamR/reference/bb_ma_plot.md)
  — all return ggplot2 objects.
- **Pipeline orchestrator**:
  [`bb_pipeline()`](https://cttir.github.io/bambamR/reference/bb_pipeline.md)
  runs the full analysis from FASTQ, BAM, or count matrix to plots in a
  single call.
- **Export**:
  [`bb_export_csv()`](https://cttir.github.io/bambamR/reference/bb_export_csv.md),
  [`bb_export_tsv()`](https://cttir.github.io/bambamR/reference/bb_export_tsv.md),
  [`bb_export_rds()`](https://cttir.github.io/bambamR/reference/bb_export_rds.md).
- **Shiny app**:
  [`bb_run_app()`](https://cttir.github.io/bambamR/reference/bb_run_app.md)
  launches an interactive analysis dashboard.
- **Example data**:
  [`bb_example_counts()`](https://cttir.github.io/bambamR/reference/bb_example_counts.md),
  [`bb_example_mutations()`](https://cttir.github.io/bambamR/reference/bb_example_mutations.md),
  and
  [`bb_example_de()`](https://cttir.github.io/bambamR/reference/bb_example_de.md)
  provide bundled datasets for exploring all features.
- Two operating modes: minimal (CRAN-only) and full (+ Bioconductor).
  All Bioconductor dependencies are optional and checked at runtime.

### Bug fixes

- None yet.
