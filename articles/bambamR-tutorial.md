# 5-Minute Tutorial: RNA-Seq Analysis with bambamR

![bambamR logo](figures/logo.svg)![bambamR logo](figures/logo-dark.svg)

[![R-CMD-check](https://github.com/CTTIR/bambamR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/CTTIR/bambamR/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/CTTIR/bambamR/actions/workflows/pkgdown.yaml/badge.svg)](https://cttir.github.io/bambamR/)
[![CRAN
status](https://www.r-pkg.org/badges/version/bambamR)](https://CRAN.R-project.org/package=bambamR)
[![Codecov test
coverage](https://codecov.io/gh/CTTIR/bambamR/branch/main/graph/badge.svg)](https://app.codecov.io/gh/CTTIR/bambamR?branch=main)
[![CRAN
downloads](https://cranlogs.r-pkg.org/badges/bambamR)](https://cran.r-project.org/package=bambamR)
[![CRAN downloads
total](https://cranlogs.r-pkg.org/badges/grand-total/bambamR)](https://cran.r-project.org/package=bambamR)
[![License:
MIT](https://img.shields.io/badge/license-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)

This tutorial walks through a complete RNA-seq analysis in under 5
minutes using only bambamR’s built-in example data. **No Bioconductor
packages are required.**

## Load the Package

``` r

library(bambamR)
#> bambamR: Optional Bioconductor engines are installed
```

## 1. Load Data and Normalize

``` r

# Bundled example: 200 genes x 10 samples (5 control, 5 treatment)
ex <- bb_example_counts()
cpm <- bb_normalize(ex$counts, method = "cpm")
cat(nrow(cpm), "genes,", ncol(cpm), "samples\n")
#> 200 genes, 10 samples
```

## 2. PCA: Do the Groups Separate?

``` r

bb_pca(cpm, ex$metadata, color_by = "condition", label = TRUE)
```

![PCA plot of the example samples coloured by condition, with each point
labelled by its sample name; control and treatment samples separate
along the first principal
component.](bambamR-tutorial_files/figure-html/pca-1.png)

## 3. Volcano Plot: What Changed?

Use the pre-computed DE results (no Bioconductor needed):

``` r

de <- bb_example_de()
bb_volcano(de, fc_cutoff = 1, p_cutoff = 0.05, n_label = 6)
```

![Volcano plot of log2 fold change against negative log10 adjusted
p-value; up- and down-regulated genes are coloured and the most
significant genes are
labelled.](bambamR-tutorial_files/figure-html/volcano-1.png)

## 4. MA Plot: Effect vs. Expression

``` r

bb_ma_plot(de)
```

![MA plot of log2 fold change against log10 mean expression with
significant genes
highlighted.](bambamR-tutorial_files/figure-html/ma-1.png)

## 5. Heatmap: Top DE Genes

``` r

bb_heatmap(cpm, de_result = de, n_genes = 20)
```

![Clustered heatmap of row-scaled expression (Z-scores) for the top
differentially expressed genes across the example
samples.](bambamR-tutorial_files/figure-html/heatmap-1.png)

## 6. Oncoplot: Mutation Landscape

``` r

mut <- bb_example_mutations()
bb_oncoplot(mut$mutations, n_genes = 10, annotation_df = mut$clinical)
```

![Oncoplot of the ten most frequently mutated genes: tiles coloured by
mutation type, bar charts of mutated genes per sample and of samples per
gene, and clinical annotation tracks
below.](bambamR-tutorial_files/figure-html/oncoplot-1.png)

## 7. Export

``` r

bb_export_csv(de, "my_results.csv")
```

## Done!

In 7 steps you went from a count matrix to publication-ready PCA,
volcano, MA, heatmap, and onco plots – all without Bioconductor.

### What’s next?

- With Bioconductor:
  [`bb_deseq2()`](https://cttir.github.io/bambamR/reference/bb_deseq2.md),
  [`bb_edger()`](https://cttir.github.io/bambamR/reference/bb_edger.md),
  [`bb_limma_voom()`](https://cttir.github.io/bambamR/reference/bb_limma_voom.md)
  for real DE analysis
- From raw data:
  [`bb_pipeline()`](https://cttir.github.io/bambamR/reference/bb_pipeline.md)
  handles FASTQ -\> alignment -\> counts -\> DE -\> plots
- Interactive:
  [`bb_run_app()`](https://cttir.github.io/bambamR/reference/bb_run_app.md)
  launches a Shiny dashboard
- Customization: every plot returns a ggplot2 object, so add
  `+ theme()`, `+ labs()`, `+ scale_*()` as needed
- See
  [`vignette("bambamR-quickstart")`](https://cttir.github.io/bambamR/articles/bambamR-quickstart.md)
  for the full walkthrough and
  [`vignette("bambamR-oncoplot")`](https://cttir.github.io/bambamR/articles/bambamR-oncoplot.md)
  for oncoplot customization
