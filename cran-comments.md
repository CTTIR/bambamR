## R CMD check results

0 errors | 0 warnings | 0 notes

## Test environments

* local: Windows 11, R 4.5.2
* GitHub Actions: ubuntu-latest (R release, devel, oldrel)
* GitHub Actions: macOS-latest (R release)
* GitHub Actions: windows-latest (R release)

## Bioconductor dependencies

Several packages in `Suggests` are from Bioconductor (DESeq2, edgeR, limma,
ShortRead, Rsamtools, GenomicAlignments, etc.). All Bioconductor dependencies
are optional. The package operates in "minimal mode" with CRAN-only
dependencies and gracefully degrades when Bioconductor packages are not
available. Every function that uses a Bioconductor package checks availability
with `requireNamespace()` first.

## Downstream dependencies

This is a new package with no downstream dependencies.
