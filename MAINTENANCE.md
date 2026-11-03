# Maintenance audit — 30 October 2026

Baseline: `820d8c96b15aa20e3e765fda3b16b2adf724375c`. Work is isolated
from the original checkout. Existing authors, examples and source data
remain the provenance of this project.

## Purpose and acceptance criteria

bambamR serves researchers importing RNA-seq counts or reads, producing
QC, normalization, differential expression and plots through R or a
local Shiny app. The immediate release criterion is preservation of
sample identity and explicit failure when an analysis cannot be
completed. No silent sample removal, alignment replacement, malformed
input acceptance or stale analysis display is acceptable. Minimal
CRAN-only use must have a runnable example; optional engines must be
distinguished from tested configurations.

## Prioritized roadmap

1.  **Essential: sample and calculation integrity.** Separate FASTQ
    outputs, preserve sample names into counts, require exact metadata
    matching, validate counts and TPM lengths, respect statistical
    designs and contrasts. Validate with independent numerical examples
    and multi-sample regression fixtures.
2.  **Essential: external command failures.** Quote arguments, separate
    SAM from diagnostics, check exit status and products, protect
    existing outputs, and request paired-fragment counts explicitly.
    Exercise real subprocesses with paths containing spaces and injected
    failures, plus supported native tools where available.
3.  **Everyday use: honest pipeline/app state.** Stop on failed
    requested DE, report plotting failures, clear old app results when
    inputs change, validate uploads and FASTQ records. Exercise a
    complete counts workflow and Shiny server state transitions,
    including failures.
4.  **Release evidence.** Synchronize function help, examples and
    migration notes; build/install/check a source package, test the
    installed first-use path and report optional backends and hosted CI
    separately. Full real-data validation and cross-platform
    external-tool qualification remain release gates where this machine
    cannot supply evidence.

Longer-term work: an explicit paired-end sample sheet, streaming
whole-library QC, strand-aware counting options, and workflow
resumption. These need a separate design and biological benchmark; they
are not inferred from filenames.

## Confirmed baseline findings

| Severity | Evidence/location | Impact/action |
|----|----|----|
| Critical | `R/pipeline.R` passes one alignment directory to every FASTQ; `R/alignment.R` uses fixed BAM names | Earlier samples overwritten; separate outputs and pass their exact paths onward |
| High | `R/alignment.R`, `R/counting.R`, BAM fallbacks ignore process status and pass raw shell arguments | False success, mixed SAM/diagnostics, unsafe path handling; checked subprocesses |
| High | `R/normalize.R` ignores names on TPM lengths and divides zero libraries | Incorrect abundance values; align by gene ID and reject invalid divisors |
| High | `R/de_analysis.R`, `.run_de()` intersect samples; edgeR ignores the pipeline design; limma indexes original design after contrasts | Dropped samples, wrong model or invalid coefficient; strict inputs and explicit model handling |
| High | `inst/app/app.R` retains old results on uploads and swallows missing DE engines | Stale results can be displayed/exported with new counts; invalidate derived state |
| Medium | FASTQ fallback truncates incomplete records; sampling reads the entire file | Silent malformed input and unnecessary memory use; bounded reads and validation |
| Medium | Internal counting ignores feature/attribute settings and forces GTF | Reject unsupported settings and let annotation format be detected |

## Baseline and verification

The first test attempt stopped before tests because `data.table` was
missing. After installing CRAN dependencies in a separate audit library,
the existing suite passed with 20 skips (mostly optional Bioconductor
engines). Baseline green CI and these mocks do not verify native
aligners or scientific correctness.

Completed milestones 1–4 locally. Recorded results:

| Check | Observed result |
|----|----|
| Regression baseline | 23 failing assertions/errors in the new integrity suite before repairs |
| Full test suite | 164 tests, 313 passing assertions, zero failures/errors/warnings; 8 skips for dependency-missing branches in the fully provisioned environment |
| `R CMD check --as-cran` | 0 errors, 0 warnings, 1 NOTE (new submission), including examples, donttest examples, vignette rebuilds, HTML help and indexed PDF manual |
| Native integration | minimap2 2.31 + samtools 1.25 + featureCounts 2.0.2; two samples retain expected counts; internal/external counts agree on the fixture; a paired fragment counts once |
| Minimal installation | Source archive installed into a separate library; README counts workflow completes without loading Bioconductor namespaces |
| Shiny | Server tests exercise normalization, settings changes, valid/invalid uploads, missing-engine failures and mismatched metadata; no browser end-to-end test is claimed |
| Documentation | 43 HTML pages built; 1,582 local references checked with no missing targets; heatmap/oncoplot annotation alignment visually inspected |
| Static checks | Actionlint and diff whitespace checks pass; standard package lint has zero findings when evaluated with the package namespace loaded; house-style findings remain report-only |
| Package contents | Repository Git pointer, citation YAML and draft publisher metadata excluded from source archive |

Environment: Apple Silicon macOS 27.0.1, R 4.6.1, data.table 1.18.6.1,
ggplot2 4.0.3, DESeq2 1.52.0, edgeR 4.10.7, limma 3.68.5, ShortRead
1.70.0, Rsamtools 2.28.0, GenomicAlignments 1.48.0, GenomicFeatures
1.64.0, txdbmaker 1.8.0 and shiny 1.14.0.

Validation commands (after installing the required check
tools/dependencies):

``` sh
Rscript --vanilla -e 'testthat::test_local()'
R CMD build .
R CMD check --as-cran bambamR_0.1.0.tar.gz
R CMD INSTALL bambamR_0.1.0.tar.gz
Rscript --vanilla tools/verify-native.R /new/empty/native-test-directory
```

The README gives the minimal prerequisites and exact first-use commands.
Native integration requires the named tools on PATH plus the internal
counting packages. Set `BAMBAMR_REQUIRE_ALL_ALIGNERS=true` to require
STAR and HISAT2 as well; the new Linux CI job does so. Without that
setting, absent optional aligners are explicitly reported as not
qualified.

The original manual-check failure was missing local TeX/makeindex,
resolved by isolated audit tools. pkgdown 2.2.1 emits upstream Pandoc
math-option deprecation warnings; the rendered site and local links
pass. Synthetic GTF import warns that no genome assembly version is
supplied; the fixture is deliberately synthetic.

## Remaining gates and next actions

Before publication, run the new revision on hosted CI, including native
STAR and HISAT2 and the supported R/OS matrix. Verify the intended
Zenodo record/DOI and publication date, then prepare the release/tag and
submission with those verified values. No push, release, deployment or
submission was performed in this audit. Study-specific biological
validity and large-data performance are not established by the synthetic
fixture. Paired-end sample-sheet automation, streaming QC, strand-aware
counting and the house-style/API migration remain deferred.

The audit logs, package archive, hashes, native fixture products and
verification receipt are retained outside the checkout in the CTTIR
public-maintenance evidence folder. The original working tree was
preserved.

Authoritative references checked for this audit:

- [R system2 argument quoting and exit
  status](https://www.stat.math.ethz.ch/R-manual/R-devel/library/base/html/system2.html)
- [Subread paired-fragment
  counting](https://subread.sourceforge.net/featureCounts.html)
- [limma reference
  manual](https://bioconductor.org/packages/release/bioc/manuals/limma/man/limma.pdf)

## Additional review findings and release provenance

- The source archive included a linked-worktree `.git` pointer and
  repository-only citation/publication metadata. Explicit build
  exclusions now keep those files out of the package while retaining
  them in Git.
- DOI `10.5281/zenodo.21889924` and the corresponding public Zenodo
  record both returned HTTP 404 on 30 October 2026; the GitHub release
  list was empty. `.zenodo.json` already labels itself a prepared,
  unpublished draft. The README no longer advertises a published DOI.
  The old DOI and date remain as comments in `CITATION.cff`, pending
  publisher verification, with all authors retained.
- Oncoplot percentages formerly used only samples with selected
  mutations. The denominator now includes every sample represented in
  the input mutation table. Annotation panels share the same grid and
  sample limits; final PNGs were visually inspected. Heatmaps now render
  requested annotation tracks and label unscaled values as expression.
- Attaching the package no longer eagerly loads every optional
  Bioconductor namespace. Engine functions still check their own runtime
  requirements.
- The shared house linter reports existing API names such as `padj` and
  `pvalue` and console output in examples. Renaming the public result
  schema would be a migration, not a safe maintenance edit. The existing
  report-only lint setting is unchanged; these findings remain
  explicitly deferred.
