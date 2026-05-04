test_that("bb_oncoplot works with simple data", {
  set.seed(42)
  mut_data <- data.frame(
    sample = sample(paste0("S", 1:10), 30, replace = TRUE),
    gene = sample(c("TP53", "KRAS", "PIK3CA", "PTEN", "APC"), 30,
                  replace = TRUE),
    mutation_type = sample(c("Missense_Mutation", "Nonsense_Mutation",
                             "Frame_Shift_Del"), 30, replace = TRUE),
    stringsAsFactors = FALSE
  )

  p <- bb_oncoplot(mut_data, n_genes = 5)
  expect_true(inherits(p, "gg") || inherits(p, "patchwork"))
})

test_that("bb_oncoplot works with MAF-format data", {
  maf <- data.frame(
    Hugo_Symbol = c("TP53", "KRAS", "PIK3CA"),
    Tumor_Sample_Barcode = c("S1", "S2", "S3"),
    Variant_Classification = c("Missense_Mutation", "Nonsense_Mutation",
                                "Splice_Site"),
    stringsAsFactors = FALSE
  )

  p <- bb_oncoplot(maf, n_genes = 3, show_barplot = FALSE)
  expect_s3_class(p, "gg")
})

test_that("bb_oncoplot handles multi-hit mutations", {
  dat <- data.frame(
    sample = c("S1", "S1", "S2"),
    gene = c("TP53", "TP53", "KRAS"),
    mutation_type = c("Missense_Mutation", "Nonsense_Mutation",
                      "Missense_Mutation"),
    stringsAsFactors = FALSE
  )

  p <- bb_oncoplot(dat, n_genes = 2, show_barplot = FALSE)
  expect_s3_class(p, "gg")
})

test_that("bb_oncoplot with annotation tracks", {
  set.seed(42)
  mut_data <- data.frame(
    sample = rep(paste0("S", 1:5), each = 3),
    gene = rep(c("TP53", "KRAS", "PTEN"), 5),
    mutation_type = sample(c("Missense_Mutation", "Nonsense_Mutation"), 15,
                           replace = TRUE),
    stringsAsFactors = FALSE
  )
  anno <- data.frame(
    stage = c("I", "II", "III", "I", "IV"),
    row.names = paste0("S", 1:5)
  )

  p <- bb_oncoplot(mut_data, n_genes = 3, annotation_df = anno,
                    show_barplot = FALSE)
  expect_true(inherits(p, "gg") || inherits(p, "patchwork"))
})

test_that("bb_oncoplot validates input", {
  expect_error(bb_oncoplot("not a data.frame"), "data.frame")
  expect_error(
    bb_oncoplot(data.frame(a = 1, b = 2, c = 3)),
    "columns"
  )
})

test_that("bb_oncoplot custom colors work", {
  dat <- data.frame(
    sample = c("S1", "S2"),
    gene = c("TP53", "KRAS"),
    mutation_type = c("Missense_Mutation", "Nonsense_Mutation"),
    stringsAsFactors = FALSE
  )
  custom <- c("Missense_Mutation" = "#FF0000", "Nonsense_Mutation" = "#0000FF")
  p <- bb_oncoplot(dat, mutation_colors = custom, show_barplot = FALSE)
  expect_s3_class(p, "gg")
})
