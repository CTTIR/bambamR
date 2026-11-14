## Submission

This is a new submission of `bambamR` (version 0.2.0).

## Test environments

* Local macOS (aarch64), R 4.6.1: full `R CMD check --as-cran` on the
  built tarball, including the PDF manual with TinyTeX and all suggested
  CRAN and Bioconductor packages. The 2026-11-14 candidate passed 364 test
  assertions. Eleven tests were skipped by their missing-package,
  platform or CRAN guards; the native minimap2/samtools regression ran.
* The native regression used minimap2 2.31-r1302 and samtools, with a
  deterministic 10 kb synthetic reference and single- and paired-end reads.
  It checked alignment coordinates, strands, pairing, unmapped records and
  BAM integrity, including paths containing spaces and apostrophes.
* Configured GitHub Actions matrix (r-lib/actions v2, `--as-cran`):
  macos-latest (release), windows-latest (release),
  ubuntu-latest (devel, release, oldrel-1). Exact-commit results must be
  verified before submission.

## R CMD check results

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE
  New submission

## Notes for CRAN

* Bioconductor packages are optional Suggests and are used after
  `requireNamespace()` checks. Examples, tests and vignettes that require
  an unavailable backend are skipped.
* External command-line tools (STAR, HISAT2, minimap2, samtools and
  featureCounts) are optional and listed in SystemRequirements. The tiny
  native regression runs only when minimap2 and samtools are available;
  it downloads no reference data. The `bb_align()` example is in
  `\dontrun{}` because it requires an installed aligner and user inputs.
  STAR, HISAT2 and standalone featureCounts were not qualified by this run.
* Examples loading optional Bioconductor backends use `\donttest{}` because
  backend loading alone can take several seconds. These examples were run
  during the full check above.
* `bb_pipeline()` writes below `tempdir()` unless the caller supplies an
  output directory. Examples and tests use temporary output paths.
* No published CRAN release or archival DOI is claimed by this candidate.
