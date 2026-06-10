# bambamR

**End-to-end RNA-seq processing from FASTQ to publication-ready plots.**

bambamR provides a streamlined toolkit for RNA-seq analysis covering
FASTQ/BAM import, quality control, alignment, read counting,
normalization, differential expression analysis, and publication-ready
visualizations including onco plots, volcano plots, heatmaps, and PCA.

## Two Operating Modes

bambamR never breaks if Bioconductor is missing. Every function that
needs an optional package checks first and either falls back gracefully
or gives a clear install instruction.

| Feature | Minimal (CRAN-only) | Full (+ Bioconductor) |
|----|:--:|:--:|
| CPM / TPM normalization | Yes | Yes |
| TMM / RLE normalization | – | Yes (edgeR / DESeq2) |
| Differential expression | – | Yes (DESeq2, edgeR, limma-voom) |
| PCA, volcano, heatmap, MA | Yes | Yes |
| **Onco plots** | **Yes** | **Yes** |
| FASTQ / BAM import | Base-R fallback | ShortRead / Rsamtools |
| Alignment wrappers | STAR / HISAT2 / minimap2 | STAR / HISAT2 / minimap2 |
| Read counting | featureCounts CLI | GenomicAlignments |
| Shiny interactive app | Yes | Yes |

## Installation

``` r

# Install from CRAN (when available)
install.packages("bambamR")

# Or install the development version from GitHub
# install.packages("pak")
pak::pak("cttir/bambamR")
```

Optionally install Bioconductor packages for full-mode features:

``` r

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "DESeq2", "edgeR", "limma",
  "ShortRead", "Rsamtools",
  "GenomicAlignments", "GenomicRanges"
))
```

## Bundled Example Datasets

bambamR ships with ready-to-use example data so you can explore every
feature without downloading external files or installing Bioconductor:

| Function | Contents |
|----|----|
| [`bb_example_counts()`](https://cttir.github.io/bambamR/reference/bb_example_counts.md) | RNA-seq count matrix (200 genes x 10 samples) with sample metadata (condition, batch, sex, age). First 20 genes have simulated DE between control and treatment. |
| [`bb_example_mutations()`](https://cttir.github.io/bambamR/reference/bb_example_mutations.md) | Mutation calls (300 mutations, 50 samples, 20 cancer genes) plus clinical annotations (Stage, Gender, Smoking) for oncoplot demos. |
| [`bb_example_de()`](https://cttir.github.io/bambamR/reference/bb_example_de.md) | Pre-computed DE results (500 genes) with log2FC, p-values, and base means – works without Bioconductor. |

A small FASTQ file (`example_reads.fastq`, 100 reads) and gene-length
table are also included under `inst/extdata/`.

## Quick Start

``` r

library(bambamR)

# ── 1. Load example data ───────────────────────────────────
ex <- bb_example_counts()
str(ex)
# List of 2
#  $ counts  : int [1:200, 1:10] ...
#  $ metadata: data.frame (10 obs. of 4 variables)

# ── 2. Normalize ───────────────────────────────────────────
cpm <- bb_normalize(ex$counts, method = "cpm")

# ── 3. PCA ─────────────────────────────────────────────────
bb_pca(cpm, ex$metadata, color_by = "condition")

# ── 4. Heatmap of top variable genes ──────────────────────
bb_heatmap(cpm, n_genes = 30)

# ── 5. Differential expression (requires DESeq2) ──────────
de <- bb_deseq2(ex$counts, ex$metadata)
bb_volcano(de)
bb_ma_plot(de)

# ... or use the pre-computed DE results (no Bioconductor!) ─
de <- bb_example_de()
bb_volcano(de, n_label = 10)
bb_ma_plot(de)

# ── 6. Heatmap of top DE genes ────────────────────────────
bb_heatmap(cpm, de_result = de, n_genes = 30)

# ── 7. Oncoplot ───────────────────────────────────────────
mut <- bb_example_mutations()
bb_oncoplot(mut$mutations, n_genes = 15, annotation_df = mut$clinical)

# ── 8. Read a FASTQ file ─────────────────────────────────
fq <- system.file("extdata", "example_reads.fastq", package = "bambamR")
reads <- bb_read_fastq(fq, n = 5)
reads

# ── 9. Export ─────────────────────────────────────────────
bb_export_csv(de, "de_results.csv")
```

## Full Pipeline (FASTQ to Results)

When you have FASTQ files, a genome index, and a GTF annotation, run the
entire pipeline in one call:

``` r

result <- bb_pipeline(
  fastq_dir    = "raw_reads/",
  genome_index = "STAR_index/",
  annotation   = "genes.gtf",
  sample_info  = sample_metadata,
  aligner      = "STAR",
  de_method    = "DESeq2",
  threads      = 8
)

# The result object bundles everything
result$counts       # count matrix
result$de_results   # DE data.frame
result$plots$pca    # ggplot object
result$plots$volcano
```

You can also enter the pipeline at any stage:

``` r

# From BAM files (skip alignment)
result <- bb_pipeline(bam_dir = "aligned/", annotation = "genes.gtf",
                      sample_info = meta)

# From a count matrix (skip alignment + counting)
result <- bb_pipeline(count_matrix = my_counts, sample_info = meta)
```

## Interactive Shiny App

``` r

bb_run_app()
```

Upload your count matrix and metadata through the browser, configure
normalization and DE parameters, run the analysis, and explore
interactive plots. Export results as CSV, RDS, or publication-quality
PDF.

## Function Reference

All exported functions use the `bb_` prefix and return standard R
objects (data.frames, matrices, ggplot objects).

### Import

- [`bb_read_fastq()`](https://cttir.github.io/bambamR/reference/bb_read_fastq.md)
  – read FASTQ files (ShortRead or base-R fallback)
- [`bb_read_bam()`](https://cttir.github.io/bambamR/reference/bb_read_bam.md)
  – read BAM files (Rsamtools or system samtools)
- [`bb_count_bam()`](https://cttir.github.io/bambamR/reference/bb_count_bam.md)
  – count reads in a BAM file

### Quality Control

- [`bb_qc()`](https://cttir.github.io/bambamR/reference/bb_qc.md) –
  compute QC metrics from FASTQ / BAM
- [`bb_qc_summary()`](https://cttir.github.io/bambamR/reference/bb_qc_summary.md)
  – tabular QC summary
- [`bb_plot_qc()`](https://cttir.github.io/bambamR/reference/bb_plot_qc.md)
  – per-position quality, GC, read-length plots

### Alignment & Counting

- [`bb_align()`](https://cttir.github.io/bambamR/reference/bb_align.md)
  – align reads with STAR, HISAT2, or minimap2
- [`bb_count_reads()`](https://cttir.github.io/bambamR/reference/bb_count_reads.md)
  – count reads per gene (GenomicAlignments or featureCounts)

### Normalization

- [`bb_normalize()`](https://cttir.github.io/bambamR/reference/bb_normalize.md)
  – CPM, TPM, TMM, or RLE normalization

### Differential Expression

- [`bb_deseq2()`](https://cttir.github.io/bambamR/reference/bb_deseq2.md)
  – DESeq2 wrapper (requires DESeq2)
- [`bb_edger()`](https://cttir.github.io/bambamR/reference/bb_edger.md)
  – edgeR quasi-likelihood wrapper (requires edgeR)
- [`bb_limma_voom()`](https://cttir.github.io/bambamR/reference/bb_limma_voom.md)
  – limma-voom wrapper (requires limma + edgeR)

### Visualization (all return ggplot objects)

- [`bb_oncoplot()`](https://cttir.github.io/bambamR/reference/bb_oncoplot.md)
  – publication-ready waterfall mutation plot
- [`bb_volcano()`](https://cttir.github.io/bambamR/reference/bb_volcano.md)
  – volcano plot with auto-labeling
- [`bb_heatmap()`](https://cttir.github.io/bambamR/reference/bb_heatmap.md)
  – clustered heatmap of top genes
- [`bb_pca()`](https://cttir.github.io/bambamR/reference/bb_pca.md) –
  PCA with metadata coloring
- [`bb_ma_plot()`](https://cttir.github.io/bambamR/reference/bb_ma_plot.md)
  – MA plot (mean expression vs. fold change)

### Pipeline & Export

- [`bb_pipeline()`](https://cttir.github.io/bambamR/reference/bb_pipeline.md)
  – run the full pipeline end-to-end
- [`bb_export_csv()`](https://cttir.github.io/bambamR/reference/bb_export_csv.md)
  /
  [`bb_export_tsv()`](https://cttir.github.io/bambamR/reference/bb_export_tsv.md)
  /
  [`bb_export_rds()`](https://cttir.github.io/bambamR/reference/bb_export_rds.md)
  – export results

### Example Data

- [`bb_example_counts()`](https://cttir.github.io/bambamR/reference/bb_example_counts.md)
  – example count matrix + metadata
- [`bb_example_mutations()`](https://cttir.github.io/bambamR/reference/bb_example_mutations.md)
  – example mutation data + clinical annotations
- [`bb_example_de()`](https://cttir.github.io/bambamR/reference/bb_example_de.md)
  – pre-computed DE results

### Shiny App

- [`bb_run_app()`](https://cttir.github.io/bambamR/reference/bb_run_app.md)
  – launch the interactive Shiny application

## Authors

- **Raban Heller** (maintainer) - SingleCellLab, BWK Ulm
- **Hanno Witte** - SingleCellLab, BWK Ulm
- **Konrad Steinestel** - SingleCellLab, BWK Ulm

## Citation

If you use bambamR in your research, please cite:

    Heller R, Witte H, Steinestel K (2026). bambamR: End-to-End RNA-Seq
    Processing from FASTQ to Publication-Ready Plots. R package version 0.1.0.
    https://github.com/cttir/bambamR

## Use of LLM tools

Portions of this package were prepared with assistance from large
language model tooling for narrowly defined, non-authorial tasks:
copyediting, prose smoothing, Markdown/LaTeX formatting, scaffolding of
boilerplate files (CI configs, build scripts), code refactoring. The
tools used were [Chat
AI](https://kisski.gwdg.de/leistungen/2-02-llm-service/), the LLM
service of KISSKI (GWDG), and a self-hosted **Mistral Small (24B,
Apache-2.0)** run locally via [Ollama](https://ollama.com/) and the
`ollamar` R package — local inference only, with no data sent to third
parties for the self-hosted model.

All scientific claims, methodological choices, analyses,
interpretations, and conclusions are the author’s own. No LLM-generated
text was incorporated without review and revision, and every reference
was verified against its DOI, arXiv ID, or ISBN.

## License

MIT
