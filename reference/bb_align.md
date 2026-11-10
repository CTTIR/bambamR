# Align Reads to Reference Genome

Wrapper for STAR, HISAT2, or minimap2 alignment tools. Checks that the
selected aligner is available on the system PATH before running.

## Usage

``` r
bb_align(
  fastq,
  genome_index,
  output_dir,
  aligner = c("STAR", "HISAT2", "minimap2"),
  threads = 4L,
  paired = FALSE,
  extra_args = NULL
)
```

## Arguments

- fastq:

  Character. Path(s) to FASTQ files. For paired-end, provide a character
  vector of length 2.

- genome_index:

  Character. Path to genome index directory (STAR) or index prefix
  (HISAT2, minimap2).

- output_dir:

  Character. New or empty output directory for BAM files.

- aligner:

  Character. One of `"STAR"`, `"HISAT2"`, `"minimap2"`.

- threads:

  Integer. Number of threads. Default `4`.

- paired:

  Logical. Paired-end mode. If `TRUE`, `fastq` must have 2 elements.
  Default `FALSE`.

- extra_args:

  Character or NULL. Additional arguments to pass to the aligner, one
  argument per element (for example `c("-x", "sr")`). HISAT2 optional
  arguments must not contain whitespace or shell metacharacters: its
  upstream wrapper reparses them through a shell. FASTQ and index paths
  containing these characters are handled with temporary symbolic links.

## Value

A list with components:

- bam:

  Path to the output BAM file

- stats:

  A data.frame with alignment statistics

- command:

  The exact command that was executed

- input_paths:

  For HISAT2, the original FASTQ and index paths before temporary
  aliasing. Temporary aliases are removed after execution.

## Examples

``` r
if (FALSE) { # \dontrun{
# Not run: needs an external aligner (STAR, HISAT2 or minimap2) on the
# PATH and a prebuilt genome index, neither of which can be shipped.
result <- bb_align("reads.fastq", "/path/to/index",
                   file.path(tempdir(), "alignment"), aligner = "HISAT2")
} # }
```
