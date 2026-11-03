test_that("HISAT2 input paths survive a second shell parse without source changes", {
  skip_on_os("windows")
  root <- withr::local_tempdir(pattern = "HISAT2 ' literal ; $ path ")
  index <- file.path(root, "large index,with punctuation")
  parts <- paste0(index, ".", seq_len(8L), ".ht2l")
  for (i in seq_along(parts)) writeLines(paste("index part", i), parts[i])
  reads <- file.path(root, "reads ' ; $ literal.fastq.gz")
  con <- gzfile(reads, "w")
  writeLines(c("@read", "ACGT", "+", "IIII"), con)
  close(con)
  originals <- c(parts, reads)
  hashes <- tools::md5sum(originals)
  staged <- NULL
  local_mocked_bindings(check_tool = function(tool) invisible(tool),
    .package = "bambamR")
  local_mocked_bindings(.run_tool = function(tool, args, stdout = TRUE) {
    if (tool == "hisat2") {
      prefix <- args[match("-x", args) + 1L]
      input <- args[match("-U", args) + 1L]
      staged <<- dirname(prefix)
      expect_true(all(file.exists(paste0(prefix, ".", seq_len(8L), ".ht2l"))))
      expect_match(input, "\\.fastq\\.gz$")
      connection <- gzfile(input)
      on.exit(close(connection), add = TRUE)
      expect_identical(readLines(connection), c("@read", "ACGT", "+", "IIII"))
      # Mirror the upstream wrapper's *second* shell parse, not system2 quoting.
      status <- system2("sh", c("-c", shQuote(paste("cat", input))), stdout = stdout)
      expect_identical(status, 0L)
    } else {
      writeLines("synthetic BAM", args[match("-o", args) + 1L])
    }
    character()
  }, .package = "bambamR")
  result <- bb_align(reads, index, file.path(root, "output"), aligner = "HISAT2")
  expect_true(file.exists(result$bam))
  expect_identical(result$input_paths, list(fastq = reads, genome_index = index))
  expect_false(dir.exists(staged))
  expect_identical(tools::md5sum(originals), hashes)
})

test_that("HISAT2 removes aliases after tool failure and refuses unsafe optional arguments", {
  skip_on_os("windows")
  root <- withr::local_tempdir(pattern = "HISAT2 failure ")
  index <- file.path(root, "index")
  for (i in seq_len(8L)) writeLines("index", paste0(index, ".", i, ".ht2"))
  reads <- file.path(root, "reads.fastq")
  writeLines(c("@read", "ACGT", "+", "IIII"), reads)
  staged <- NULL
  local_mocked_bindings(check_tool = function(tool) invisible(tool),
    .package = "bambamR")
  local_mocked_bindings(.run_tool = function(tool, args, stdout = TRUE) {
    staged <<- dirname(args[match("-x", args) + 1L])
    stop("synthetic aligner failure")
  }, .package = "bambamR")
  expect_error(bb_align(reads, index, file.path(root, "failed"), aligner = "HISAT2"),
    "synthetic aligner failure")
  expect_false(dir.exists(staged))
  expect_true(file.exists(reads))
  for (value in c("sample name", "; touch sentinel", "$(touch sentinel)", "`touch sentinel`")) {
    expect_error(bb_align(reads, index, tempfile(tmpdir = root), aligner = "HISAT2",
      extra_args = c("--rg-id", value)), "optional arguments cannot safely")
  }
})


test_that("Windows drive paths use safe separators without index staging", {
  helper <- .hisat2_inputs
  lookup <- new.env(parent = environment(helper))
  lookup$.Platform <- list(OS.type = "windows")
  environment(helper) <- lookup
  reads <- c("C:\\reads\\mate1.fastq", "C:\\reads\\mate2.fastq")
  prefix <- "D:\\indexes\\hg38"
  x <- helper(reads, prefix)
  expect_identical(x$fastq, c("C:/reads/mate1.fastq", "C:/reads/mate2.fastq"))
  expect_identical(x$genome_index, "D:/indexes/hg38")
  expect_null(x$directory)
  expect_identical(reads, c("C:\\reads\\mate1.fastq", "C:\\reads\\mate2.fastq"))
})

test_that("Unix literal backslashes remain filename bytes and are staged", {
  skip_on_os("windows")
  root <- withr::local_tempdir()
  reads <- file.path(root, "literal\\backslash.fastq")
  writeLines(c("@read", "ACGT", "+", "IIII"), reads)
  prefix <- file.path(root, "index")
  parts <- paste0(prefix, ".", seq_len(8L), ".ht2")
  for (part in parts) writeLines("index", part)
  before <- tools::md5sum(c(reads, parts))
  x <- .hisat2_inputs(reads, prefix)
  withr::defer(unlink(x$directory, recursive = TRUE))
  expect_true(dir.exists(x$directory))
  expect_identical(readLines(x$fastq), readLines(reads))
  expect_identical(Sys.readlink(x$fastq), normalizePath(reads))
  expect_identical(tools::md5sum(c(reads, parts)), before)
})
