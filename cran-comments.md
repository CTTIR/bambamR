## Submission

This is a new submission of `bambamR` (version 0.2.0).

## Test environments

* local macOS (aarch64), R 4.6.1:
  `R CMD check --as-cran --no-manual` (LaTeX is not installed locally, so
  the PDF manual was not built here), run twice:
  * with all suggested CRAN and Bioconductor packages installed, and
  * in a library that contains CRAN packages only
    (`_R_CHECK_FORCE_SUGGESTS_=false`), i.e. without any Bioconductor
    package.
* GitHub Actions (r-lib/actions v2, `--as-cran`):
  macos-latest (release), windows-latest (release),
  ubuntu-latest (devel, release, oldrel-1)

## R CMD check results

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE
  New submission

## Notes for CRAN

* All Bioconductor packages are in Suggests ('ShortRead', 'Biostrings',
  'Rsamtools', 'GenomicAlignments', 'GenomicFeatures', 'txdbmaker',
  'GenomicRanges', 'DESeq2', 'edgeR', 'limma', 'SummarizedExperiment',
  'Biobase', 'ComplexHeatmap'). They are only used after
  `requireNamespace()` checks; examples, tests and vignettes that need them
  are skipped when they are not installed. The package passes
  `R CMD check` with CRAN-only dependencies ("minimal mode").
* External command-line tools (STAR, HISAT2, minimap2, samtools,
  featureCounts) are optional and listed in SystemRequirements. They are
  never called during examples, tests or vignettes on CRAN; the example of
  `bb_align()` is in `\dontrun{}` because it needs an aligner and a genome
  index that cannot be shipped.
* Examples that load the optional 'Bioconductor' back ends are wrapped in
  `\donttest{}` because loading those packages alone can take several
  seconds.
* `bb_pipeline()` writes to a new directory below `tempdir()` unless the
  user supplies `output_dir`; examples, tests and vignettes only write to
  `tempdir()`.
* README.md and the vignettes carry CRAN badges linking to
  <https://CRAN.R-project.org/package=bambamR>, which returns 404 until the
  package has been accepted; all other URLs resolve.
