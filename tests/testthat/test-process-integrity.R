test_that("subprocess arguments remain literal and diagnostics stay separate", {
  dir <- withr::local_tempdir(pattern = "space path ")
  script <- file.path(dir, "echo args.R")
  writeLines(c('message("diagnostic only")',
               'cat(commandArgs(trailingOnly = TRUE), sep = "\\n")'), script)
  literal <- c("two words", "semi;colon", "$(echo unsafe)", "a'quote", "back`tick")
  result <- bambamR:::.run_tool(file.path(R.home("bin"), "Rscript"),
                                c("--vanilla", script, literal))
  expect_identical(result, literal)
  target <- file.path(dir, "stdout.txt")
  bambamR:::.run_tool(file.path(R.home("bin"), "Rscript"),
                      c("--vanilla", script, "SAM"), stdout = target)
  expect_equal(readLines(target), "SAM")
  writeLines('message("invalid index"); quit(status = 23L)', script)
  expect_error(bambamR:::.run_tool(file.path(R.home("bin"), "Rscript"),
                                   c("--vanilla", script)), "23.*invalid index")
})

test_that("alignment refuses existing products and missing successful output", {
  dir <- withr::local_tempdir()
  fq <- withr::local_tempfile()
  writeLines(c("@r", "AC", "+", "II"), fq)
  existing <- file.path(dir, "existing.bam")
  writeLines("preserve", existing)
  expect_error(bb_align(fq, "index", dir), "must be empty")
  expect_equal(readLines(existing), "preserve")
  local_mocked_bindings(check_tool = function(tool) tool,
                        .run_tool = function(...) character(), .package = "bambamR")
  expect_error(bambamR:::.align_star(fq, "index", withr::local_tempdir(),
                                     1L, FALSE, NULL), "non-empty output")
})

test_that("featureCounts requests fragments and rejects failed or malformed output", {
  args_used <- NULL
  write_counts <- function(tool, args, ...) {
    args_used <<- args
    output <- args[match("-o", args) + 1L]
    writeLines(c("# featureCounts", "Geneid\tChr\tStart\tEnd\tStrand\tLength\ta.bam",
                 "g1\tchr1\t1\t5\t+\t5\t2"), output)
    character()
  }
  local_mocked_bindings(check_tool = identity, .run_tool = write_counts,
                        .package = "bambamR")
  counts <- bambamR:::.count_featurecounts("a.bam", "ann.gtf", 1L, "exon", "gene_id", TRUE)
  expect_equal(unname(counts[1, 1]), 2)
  expect_true(all(c("-p", "--countReadPairs") %in% args_used))
  expect_error(bambamR:::.count_featurecounts(c("a.bam", "b.bam"), "ann.gtf", 1L,
                                              "exon", "gene_id", FALSE), "Invalid featureCounts")
})

test_that("named BAM paths preserve explicit sample identities", {
  expect_equal(bambamR:::.clean_sample_names(c(a = "a/aligned.bam", b = "b/aligned.bam")),
               c("a", "b"))
  expect_error(bambamR:::.clean_sample_names(c("a/aligned.bam", "b/aligned.bam")), "unique")
})
