# bambamR 0.1.0

## New features

* Initial release of bambamR.
* **Data import**: `bb_read_fastq()` and `bb_read_bam()` with ShortRead/Rsamtools
  or base-R/system-samtools fallback.
* **Quality control**: `bb_qc()`, `bb_qc_summary()`, and `bb_plot_qc()` for
  per-read quality, GC content, and read-length distributions.
* **Alignment**: `bb_align()` wrapping STAR, HISAT2, and minimap2.
* **Read counting**: `bb_count_reads()` via GenomicAlignments or featureCounts.
* **Normalization**: `bb_normalize()` supporting CPM, TPM (base R), TMM (edgeR),
  and RLE (DESeq2).
* **Differential expression**: `bb_deseq2()`, `bb_edger()`, and
  `bb_limma_voom()` returning standardized result data.frames.
* **Visualization**: `bb_oncoplot()`, `bb_volcano()`, `bb_heatmap()`,
  `bb_pca()`, `bb_ma_plot()` — all return ggplot2 objects.
* **Pipeline orchestrator**: `bb_pipeline()` runs the full analysis from FASTQ,
  BAM, or count matrix to plots in a single call.
* **Export**: `bb_export_csv()`, `bb_export_tsv()`, `bb_export_rds()`.
* **Shiny app**: `bb_run_app()` launches an interactive analysis dashboard.
* **Example data**: `bb_example_counts()`, `bb_example_mutations()`, and
  `bb_example_de()` provide bundled datasets for exploring all features.
* Two operating modes: minimal (CRAN-only) and full (+ Bioconductor).
  All Bioconductor dependencies are optional and checked at runtime.

## Bug fixes

* None yet.
