# Tests for the external-tool shell-out paths in R/alignment.R, R/counting.R,
# and R/import_bam.R. These code paths normally invoke STAR/HISAT2/minimap2,
# featureCounts, and samtools. We mock the package-internal check_tool() guard
# (so the tool "exists") and base::system2 (so nothing is actually executed),
# which lets us exercise command construction and output parsing deterministically.

test_that(".read_bam_system parses tab-delimited SAM output", {
  sam_lines <- c(
    paste("r1", "0", "chr1", "100", "60", "50M",
          "*", "0", "0", "ACGT", "IIII", sep = "\t"),
    paste("r2", "16", "chr2", "250", "30", "40M",
          "*", "0", "0", "TGCA", "HHHH", sep = "\t")
  )
  local_mocked_bindings(check_tool = function(tool) invisible("samtools"),
                        .package = "bambamR")
  local_mocked_bindings(system2 = function(...) sam_lines, .package = "base")

  df <- bambamR:::.read_bam_system("dummy.bam")
  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 2L)
  expect_equal(df$qname, c("r1", "r2"))
  expect_equal(df$flag, c(0L, 16L))
  expect_equal(df$rname, c("chr1", "chr2"))
  expect_equal(df$pos, c(100L, 250L))
})

test_that(".read_bam_system returns an empty frame for empty output", {
  local_mocked_bindings(check_tool = function(tool) invisible("samtools"),
                        .package = "bambamR")
  local_mocked_bindings(system2 = function(...) character(0), .package = "base")
  df <- bambamR:::.read_bam_system("dummy.bam")
  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 0L)
})

test_that("bb_count_bam uses system samtools when Rsamtools is unavailable", {
  local_mocked_bindings(requireNamespace = function(pkg, ...) {
    if (identical(pkg, "Rsamtools")) FALSE else TRUE
  }, .package = "base")
  local_mocked_bindings(check_tool = function(tool) invisible("samtools"),
                        .package = "bambamR")
  local_mocked_bindings(system2 = function(...) "  42  ", .package = "base")

  bam <- withr::local_tempfile(fileext = ".bam")
  writeLines("x", bam)
  expect_equal(bb_count_bam(bam), 42L)
})

test_that("bb_read_bam routes to the system fallback without Rsamtools", {
  local_mocked_bindings(requireNamespace = function(pkg, ...) {
    if (identical(pkg, "Rsamtools")) FALSE else TRUE
  }, .package = "base")
  local_mocked_bindings(check_tool = function(tool) invisible("samtools"),
                        .package = "bambamR")
  local_mocked_bindings(
    system2 = function(...) paste("r1", "0", "chr1", "100", "60", "50M",
                                  "*", "0", "0", "ACGT", "IIII", sep = "\t"),
    .package = "base"
  )
  bam <- withr::local_tempfile(fileext = ".bam")
  writeLines("x", bam)
  expect_message(df <- bb_read_bam(bam), "system samtools")
  expect_equal(nrow(df), 1L)
})

test_that(".align_star builds a command and parses the STAR log", {
  fq <- withr::local_tempfile(fileext = ".fastq")
  writeLines(c("@r1", "ACGT", "+", "IIII"), fq)
  out_dir <- withr::local_tempdir()

  # Pre-create the STAR Log.final.out the parser will read.
  prefix <- file.path(out_dir, "STAR_")
  writeLines(c("Number of input reads |\t100"),
             paste0(prefix, "Log.final.out"))

  local_mocked_bindings(check_tool = function(tool) invisible(tool),
                        .package = "bambamR")
  local_mocked_bindings(system2 = function(...) "ok", .package = "base")

  res <- bambamR:::.align_star(fq, "idx", out_dir, threads = 2L,
                               paired = FALSE, extra_args = NULL)
  expect_type(res, "list")
  expect_true(all(c("bam", "stats", "command") %in% names(res)))
  expect_match(res$command, "^STAR ")
  expect_match(res$bam, "Aligned.sortedByCoord.out.bam$")
})

test_that(".align_star adds zcat for gzipped input and extra_args", {
  fq <- withr::local_tempfile(fileext = ".fastq.gz")
  con <- gzfile(fq, "w"); writeLines(c("@r1", "ACGT", "+", "IIII"), con)
  close(con)
  out_dir <- withr::local_tempdir()

  local_mocked_bindings(check_tool = function(tool) invisible(tool),
                        .package = "bambamR")
  local_mocked_bindings(system2 = function(...) "ok", .package = "base")

  res <- bambamR:::.align_star(fq, "idx", out_dir, threads = 1L,
                               paired = FALSE, extra_args = c("--quantMode",
                                                              "GeneCounts"))
  expect_match(res$command, "zcat")
  expect_match(res$command, "GeneCounts")
})

test_that(".align_hisat2 builds a paired-end command", {
  fq1 <- withr::local_tempfile(fileext = "_1.fastq")
  fq2 <- withr::local_tempfile(fileext = "_2.fastq")
  writeLines(c("@r1", "ACGT", "+", "IIII"), fq1)
  writeLines(c("@r1", "ACGT", "+", "IIII"), fq2)
  out_dir <- withr::local_tempdir()

  local_mocked_bindings(check_tool = function(tool) invisible(tool),
                        .package = "bambamR")
  local_mocked_bindings(system2 = function(...) "ok", .package = "base")

  res <- bambamR:::.align_hisat2(c(fq1, fq2), "idx", out_dir, threads = 2L,
                                 paired = TRUE, extra_args = NULL)
  expect_match(res$command, "hisat2")
  expect_match(res$command, "samtools sort")
  expect_equal(res$stats$metric, "alignment_complete")
})

test_that(".align_minimap2 builds a command and reports completion", {
  fq <- withr::local_tempfile(fileext = ".fastq")
  writeLines(c("@r1", "ACGT", "+", "IIII"), fq)
  out_dir <- withr::local_tempdir()

  local_mocked_bindings(check_tool = function(tool) invisible(tool),
                        .package = "bambamR")
  local_mocked_bindings(system2 = function(...) invisible(0L), .package = "base")

  res <- bambamR:::.align_minimap2(fq, "idx", out_dir, threads = 4L,
                                   paired = FALSE, extra_args = "-x sr")
  expect_match(res$command, "minimap2")
  expect_match(res$bam, "minimap2_aligned.bam$")
})

test_that(".count_featurecounts parses a featureCounts table", {
  out_dir <- withr::local_tempdir()
  bam1 <- file.path(out_dir, "sampleA.bam")
  bam2 <- file.path(out_dir, "sampleB.bam")
  file.create(bam1, bam2)

  # featureCounts writes a header line, a column-name line, then data.
  fc_table <- rbind(
    c("Geneid", "Chr", "Start", "End", "Strand", "Length", bam1, bam2),
    c("G1", "chr1", "100", "200", "+", "101", "10", "20"),
    c("G2", "chr1", "300", "400", "+", "101", "5", "15")
  )

  local_mocked_bindings(check_tool = function(tool) invisible(tool),
                        .package = "bambamR")
  local_mocked_bindings(
    system2 = function(command, args, ...) {
      # locate the -o output path and write the featureCounts table there
      o_idx <- which(args == "-o")
      out_file <- args[o_idx + 1L]
      writeLines("# Program:featureCounts", out_file)
      con <- file(out_file, open = "a")
      utils::write.table(fc_table, con, sep = "\t", quote = FALSE,
                         row.names = FALSE, col.names = FALSE)
      close(con)
      invisible(0L)
    },
    .package = "base"
  )

  counts <- bambamR:::.count_featurecounts(c(bam1, bam2), "ann.gtf",
                                           threads = 2L, feature_type = "exon",
                                           attr_type = "gene_id", paired = FALSE)
  expect_true(is.matrix(counts))
  expect_equal(rownames(counts), c("G1", "G2"))
  expect_equal(colnames(counts), c("sampleA", "sampleB"))
  expect_equal(unname(counts["G1", ]), c(10, 20))
})
