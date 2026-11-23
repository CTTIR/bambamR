# bambamR 0.2.0

## Bug fixes

* `bb_align()` with `aligner = "STAR"` now decompresses gzipped FASTQ files
  with `gunzip -c` instead of `zcat`. The `zcat` shipped with macOS only reads
  `.Z` files, so STAR could not read `.fastq.gz` input there.
* `bb_edger()` now stops with a clear message when `group` has a single level
  and no `design` is given, or when `group` contains missing values. Both
  cases previously failed with an unrelated error from the model matrix.
* `bb_qc()` computes GC content and per-position quality in time linear in the
  number of bases; the per-position summary previously looped over every read
  for every position. Reads are processed in blocks of about ten million
  bases, so memory stays bounded for long reads as well. Results are
  unchanged.
* `bb_qc()`: `use_fastqc` is documented as reserved and ignored; it never
  called FastQC.
* `bb_oncoplot()` no longer fails when a custom `mutation_colors` vector has no
  `"Other"` entry and the data contain an unlisted mutation type (a neutral
  grey is used), or when `genes` contains duplicates. `n_genes`, `genes`,
  `mutation_colors` and `annotation_df` are validated, and a warning is raised
  when no row name of `annotation_df` matches a sample instead of dropping the
  tracks silently.
* `bb_edger()`, `bb_limma_voom()` and TMM normalization no longer print the
  edgeR 4.0 notice that `calcNormFactors()` was renamed to `normLibSizes()`.
* `bb_volcano()` no longer warns about removed rows for every unlabelled gene,
  and drops rows without a fold change. `bb_volcano()` and `bb_ma_plot()`
  accept an unnamed `colors` vector in the documented order and report wrong
  names clearly; `n_label = 0` draws no labels.
* `bb_heatmap()` shows a single-sample matrix as zeros instead of missing
  tiles when scaling, and `bb_heatmap()` and `bb_pca()` validate `n_genes`
  and the metadata arguments.

## Plot readability

* `bb_oncoplot()`: the legend is collected below the whole figure with its
  title on top and at most three columns, so it is no longer clipped at the
  device edge; the top barplot title is split over two lines and no longer
  runs into the plot title, and its axis of gene counts uses whole-number
  breaks.
  Annotation tracks drawn without barplots collect their legends in the same
  way. When moving the legend, reset its layout too: for the patchwork
  composition use `p & ggplot2::theme(legend.position = "right") &
  ggplot2::guides(fill = ggplot2::guide_legend(ncol = 1))`, for a single
  ggplot (`show_barplot = FALSE`, no `annotation_df`) use `+` instead of `&`.
* `bb_heatmap()`: the Z-score colour scale is symmetric, so the middle colour
  always marks zero; with `annotation_col`, sample labels are drawn once below
  the annotation track instead of twice and overlapping.
* `bb_plot_qc()`: the Q20 and Q30 thresholds are labelled and use a
  colour-blind-safe pair with different line types; the read-length histogram
  uses one bin per length for narrow ranges and whole-number axis breaks. The
  combined panel draws the file legend once instead of three times, which
  leaves more room for the panels.

* Shiny app: plots are drawn at 96 dpi instead of 72 dpi, so axis and gene
  labels are readable in the result tabs. The PCA and MA tabs say why they are
  empty (no `condition` column, or a DE method without mean expression)
  instead of staying blank, and the DE table shows three significant digits.

## Documentation

* Documentation uses British English (`Language: en-GB`); code identifiers,
  argument names and plot labels are unchanged.
* The README citation, `CITATION.cff` and `.zenodo.json` name the author as
  R. Heller and use the package title of `citation("bambamR")`.
* Figures in the tutorial, quick-start and oncoplot vignettes carry alt text.
* The custom volcano example in the vignettes uses cut-off values that leave
  significant genes, so the custom colours are visible.
* The GitHub icon in the site navigation bar has an accessible label.

## CRAN preparation

* `bb_pipeline()` now writes to a new directory below `tempdir()` by default
  instead of `"bambamR_output"` in the working directory; pass `output_dir`
  to keep results.
* Examples no longer contain commented-out calls: BAM and counting examples
  build a tiny BAM file when 'Rsamtools' is installed, the pipeline example
  runs in minimal mode, and `bb_align()` documents why its example cannot run.
* Documented return values of the print and summary methods.
* `DESCRIPTION`: shorter title, quoted software names, `Language` and
  `SystemRequirements` for the optional external tools.

## Integrity fixes

* Preserve requested BAM columns and their types when the Rsamtools reader
  encounters a valid file with no alignments. Empty results now support the
  same field selection as populated results.

* Expand existing Windows short read-path aliases before HISAT2 admission;
  preserve bare index names and the caller's original path provenance.

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
argument tokens. limma returns one contrast per call.

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
