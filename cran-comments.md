## Local R CMD check results (30 October 2026)

0 errors | 0 warnings | 1 note

The note is CRAN incoming feasibility: "New submission".
The full check includes examples, donttest examples, tests, vignette rebuilding,
HTML help validation and the indexed PDF manual.

## Verified environment

* macOS 27.0.1 (Apple Silicon), R 4.6.1.
* CRAN-only first use from a separate package library.
* Optional Bioconductor engines, including DESeq2, edgeR, limma,
  ShortRead, Rsamtools and GenomicAlignments.
* Native synthetic-read integration with minimap2 2.31, samtools 1.25 and
  featureCounts 2.0.2; counts checked against known expected values.

The new revision has not run on hosted GitHub Actions, Windows or Linux.
Native STAR/HISAT2 integration is prepared in CI but was not run on this Mac.
These must be qualified before claiming the complete supported-platform matrix.

## Publication status

No CRAN submission or GitHub release was performed. The previously recorded
Zenodo identifier does not currently resolve; its provenance is preserved in
CITATION.cff comments and MAINTENANCE.md. Confirm publisher metadata before release.

## Optional dependencies

Statistical engines and read-processing backends are optional. Requested missing
engines fail with installation guidance. CRAN-only pipeline use explicitly skips
DE; it does not report a failed DE analysis as complete.

See MAINTENANCE.md for the audit, migration notes, exact local validation commands,
and remaining qualification work.
