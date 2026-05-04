test_that("bb_export_rds round-trips correctly", {
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp))

  dat <- data.frame(gene = c("A", "B"), log2fc = c(1, -1))
  bb_export_rds(dat, tmp)
  expect_true(file.exists(tmp))
  expect_equal(readRDS(tmp), dat)
})

test_that("bb_export_csv writes correct format", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))

  de <- data.frame(
    gene = paste0("gene", 1:3),
    log2fc = c(1.5, -2.0, 0.3),
    pvalue = c(0.001, 0.01, 0.5),
    padj = c(0.01, 0.05, 0.8)
  )
  bb_export_csv(de, tmp)
  expect_true(file.exists(tmp))

  reimport <- data.table::fread(tmp, data.table = FALSE)
  expect_equal(reimport$gene, de$gene)
  expect_equal(reimport$log2fc, de$log2fc, tolerance = 1e-6)
})

test_that("bb_export_tsv writes tab-delimited", {
  tmp <- tempfile(fileext = ".tsv")
  on.exit(unlink(tmp))

  de <- data.frame(gene = "A", log2fc = 1, pvalue = 0.01, padj = 0.05)
  bb_export_tsv(de, tmp)
  expect_true(file.exists(tmp))

  lines <- readLines(tmp)
  expect_true(grepl("\t", lines[2]))
})

test_that("export validates path argument", {
  expect_error(bb_export_rds(1), "path")
  expect_error(bb_export_csv(data.frame(a = 1)), "path")
  expect_error(bb_export_tsv(data.frame(a = 1)), "path")
})
