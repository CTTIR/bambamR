test_that("native minimap2 and samtools preserve synthetic alignments", {
  skip_if(!nzchar(Sys.which("minimap2")), "minimap2 is not installed")
  skip_if(!nzchar(Sys.which("samtools")), "samtools is not installed")
  parent <- withr::local_tempdir()
  work <- file.path(parent, "fixtures' space")
  dir.create(work)
  withr::local_seed(123)
  reference <- paste0(sample(c("A", "C", "G", "T"), 10000L,
                             replace = TRUE), collapse = "")
  fasta <- file.path(work, "reference's genome.fa")
  writeLines(c(">chrSynthetic", reference), fasta)
  index <- file.path(work, "reference's index.mmi")
  bambamR:::.run_tool("minimap2", c("-x", "sr", "-d", index, fasta))
  revcomp <- function(x) {
    paste0(rev(strsplit(chartr("ACGT", "TGCA", x), "")[[1L]]), collapse = "")
  }
  write_fastq <- function(path, ids, sequences) {
    records <- Map(function(id, sequence) {
      c(paste0("@", id), sequence, "+", strrep("I", nchar(sequence)))
    }, ids, sequences)
    writeLines(unlist(records, use.names = FALSE), path)
    path
  }
  fastq <- write_fastq(file.path(work, "single's reads.fastq"),
                      c("forward", "reverse", "unmapped"),
                      c(substr(reference, 1001L, 1150L),
                        revcomp(substr(reference, 2501L, 2650L)),
                        strrep("N", 150L)))
  single <- suppressMessages(bb_align(
    fastq, index, file.path(work, "single's output"),
    aligner = "minimap2", threads = 1L, extra_args = c("-x", "sr")
  ))
  expect_silent(bambamR:::.run_tool("samtools", c("quickcheck", single$bam)))
  reads <- bambamR:::.read_bam_system(single$bam)
  expect_equal(nrow(reads), 3L)
  expect_setequal(reads$qname, c("forward", "reverse", "unmapped"))
  reads <- reads[match(c("forward", "reverse", "unmapped"), reads$qname), ]
  expect_identical(reads$pos, c(1001L, 2501L, 0L))
  expect_identical(reads$cigar, c("150M", "150M", "*"))
  expect_equal(bitwAnd(reads$flag, c(4L, 16L, 4L)), c(0L, 16L, 4L))

  first <- write_fastq(file.path(work, "mate one.fastq"), "pair1",
                       substr(reference, 1001L, 1150L))
  second <- write_fastq(file.path(work, "mate two.fastq"), "pair1",
                        revcomp(substr(reference, 1251L, 1400L)))
  paired <- suppressMessages(bb_align(
    c(first, second), index, file.path(work, "paired's output"),
    aligner = "minimap2", threads = 1L, paired = TRUE,
    extra_args = c("-x", "sr")
  ))
  expect_silent(bambamR:::.run_tool("samtools", c("quickcheck", paired$bam)))
  pairs <- bambamR:::.read_bam_system(paired$bam)
  expect_equal(nrow(pairs), 2L)
  expect_identical(pairs$pos, c(1001L, 1251L))
  expect_true(all(bitwAnd(pairs$flag, 2L) == 2L))
  expect_identical(pairs$cigar, c("150M", "150M"))
})
