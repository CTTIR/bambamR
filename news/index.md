# Changelog

## bambamR 0.2.0

### CRAN preparation

- [`bb_pipeline()`](https://cttir.github.io/bambamR/reference/bb_pipeline.md)
  now writes to a new directory below
  [`tempdir()`](https://rdrr.io/r/base/tempfile.html) by default instead
  of `"bambamR_output"` in the working directory; pass `output_dir` to
  keep results.
- Examples no longer contain commented-out calls: BAM and counting
  examples build a tiny BAM file when ‘Rsamtools’ is installed, the
  pipeline example runs in minimal mode, and
  [`bb_align()`](https://cttir.github.io/bambamR/reference/bb_align.md)
  documents why its example cannot run.
- Documented return values of the print and summary methods.
- `DESCRIPTION`: shorter title, quoted software names, `Language` and
  `SystemRequirements` for the optional external tools.

### Integrity fixes

- Expand existing Windows short read-path aliases before HISAT2
  admission; preserve bare index names and the caller’s original path
  provenance.

- Normalize Windows path separators before HISAT2 alias admission while
  preserving literal backslashes in Unix filenames and original input
  paths.

- Preserve HISAT2 FASTQ/index paths containing spaces or shell
  punctuation with private temporary aliases, including compressed reads
  and large indexes. Retain original input paths in the result; reject
  unsafe optional arguments that the upstream wrapper would parse
  through a second shell.

- Preserve separate FASTQ alignment products and sample names; reject
  ambiguous pipeline entry points, incompatible skipped steps and output
  replacement.

- Check external command status/products and quote argument tokens.
  HISAT2 SAM output no longer contains diagnostics. Paired featureCounts
  requests fragments.

- Match named TPM lengths and metadata by identifier, reject invalid
  counts and empty libraries, preserve edgeR adjustment variables, and
  repair limma contrasts.

- Fail requested DE explicitly; retain plot warnings in the saved
  result.

- Validate FASTQ records and limit reads when a sample size is
  requested.

- Clear stale Shiny results, implement previews and normalization
  selection, and provide a normalization-only mode. TPM remains
  available through the R API.

- Preserve the oncoplot cohort denominator, display requested annotation
  tracks, and label unscaled heatmaps as expression rather than
  Z-scores.

Compatibility: use a fresh output directory for each alignment/pipeline
run. Metadata must now match all samples. Supply external-tool options
as separate argument tokens. limma returns one contrast per call. See
MAINTENANCE.md for validation evidence and remaining release gates.

## bambamR 0.1.0

- Standardize the named v5 package logo across README, documentation
  pages and available application headers, with matching light and dark
  artwork.

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
