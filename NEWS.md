# Development changes

* Normalize Windows path separators before HISAT2 alias admission while
  preserving literal backslashes in Unix filenames and original input paths.
* Preserve HISAT2 FASTQ/index paths containing spaces or shell punctuation with
  private temporary aliases, including compressed reads and large indexes.
  Retain original input paths in the result; reject unsafe optional arguments
  that the upstream wrapper would parse through a second shell.

* Preserve separate FASTQ alignment products and sample names; reject ambiguous
  pipeline entry points, incompatible skipped steps and output replacement.
* Check external command status/products and quote argument tokens. HISAT2 SAM
  output no longer contains diagnostics. Paired featureCounts requests fragments.
* Match named TPM lengths and metadata by identifier, reject invalid counts and
  empty libraries, preserve edgeR adjustment variables, and repair limma contrasts.
* Fail requested DE explicitly; retain plot warnings in the saved result.
* Validate FASTQ records and limit reads when a sample size is requested.
* Clear stale Shiny results, implement previews and normalization selection, and
  provide a normalization-only mode. TPM remains available through the R API.
* Preserve the oncoplot cohort denominator, display requested annotation tracks,
  and label unscaled heatmaps as expression rather than Z-scores.

Compatibility: use a fresh output directory for each alignment/pipeline run.
Metadata must now match all samples. Supply external-tool options as separate
argument tokens. limma returns one contrast per call. See MAINTENANCE.md for
validation evidence and remaining release gates.

# bambamR 0.1.0

* Standardize the named v5 package logo across README, documentation pages
  and available application headers, with matching light and dark artwork.


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
