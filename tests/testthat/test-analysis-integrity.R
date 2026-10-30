integrity_counts <- function() {
  matrix(c(10, 20, 30, 40), 2,
         dimnames = list(c("g1", "g2"), c("s1", "s2")))
}

test_that("TPM matches named lengths by gene rather than vector position", {
  counts <- integrity_counts()
  expected <- bb_normalize(counts, "tpm", c(1000, 2000))
  expect_equal(bb_normalize(counts, "tpm", c(g2 = 2000, g1 = 1000)), expected)
  expect_error(bb_normalize(counts, "tpm", c(g1 = 1000, other = 2000)), "gene")
  for (lengths in list(c(0, 1000), c(-1, 1000), c(NA, 1000), c(Inf, 1000))) {
    expect_error(bb_normalize(counts, "tpm", lengths), "positive|finite")
  }
})

test_that("invalid counts and zero libraries cannot produce plausible output", {
  for (value in c(-1, NA, Inf)) {
    counts <- integrity_counts()
    counts[1, 1] <- value
    expect_error(bb_normalize(counts), "finite|negative")
  }
  counts <- integrity_counts()
  counts[, 1] <- 0
  expect_error(bb_normalize(counts), "library")
  counts <- integrity_counts()
  colnames(counts) <- c("duplicate", "duplicate")
  expect_error(bb_normalize(counts), "unique")
})

test_that("metadata mismatches cannot silently remove samples", {
  counts <- integrity_counts()
  meta <- data.frame(condition = c("a", "b"), row.names = c("s1", "other"))
  expect_error(bb_pipeline(count_matrix = counts, sample_info = meta,
                           output_dir = withr::local_tempdir(), skip = c("de", "viz")),
               "match")
})

test_that("pipeline preserves separate alignment products and sample IDs", {
  dir <- withr::local_tempdir()
  fq_dir <- file.path(dir, "reads")
  dir.create(fq_dir)
  for (name in c("s1", "s2")) {
    writeLines(c(paste0("@", name), "AC", "+", "II"),
               file.path(fq_dir, paste0(name, ".fastq")))
  }
  paths <- character()
  local_mocked_bindings(
    bb_align = function(fastq, genome_index, output_dir, ...) {
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
      bam <- file.path(output_dir, "aligned.bam")
      writeLines(basename(fastq), bam)
      paths <<- c(paths, bam)
      list(bam = bam, stats = data.frame(metric = "reads", value = "1"))
    },
    bb_count_reads = function(bam_paths, ...) {
      expect_length(unique(bam_paths), 2)
      expect_equal(vapply(bam_paths, readLines, character(1), USE.NAMES = FALSE),
                   c("s1.fastq", "s2.fastq"))
      counts <- integrity_counts()
      colnames(counts) <- basename(bam_paths)
      counts
    }, .package = "bambamR"
  )
  meta <- data.frame(condition = c("a", "b"), row.names = c("s1", "s2"))
  result <- bb_pipeline(fastq_dir = fq_dir, genome_index = "index",
                        annotation = "genes.gtf", sample_info = meta,
                        output_dir = file.path(dir, "output"),
                        skip = c("qc", "de", "viz"))
  expect_length(unique(paths), 2)
  expect_equal(colnames(result$counts), rownames(meta))
})

test_that("failed requested DE never reports a complete pipeline", {
  local_mocked_bindings(.run_de = function(...) stop("model failed"),
                        .package = "bambamR")
  meta <- data.frame(condition = c("a", "b"), row.names = c("s1", "s2"))
  out <- withr::local_tempdir()
  expect_error(bb_pipeline(count_matrix = integrity_counts(), sample_info = meta,
                           output_dir = out, skip = "viz"), "model failed")
  expect_false(file.exists(file.path(out, "bambamR_result.rds")))
})

test_that("a failed STAR process does not return a BAM path as success", {
  local_mocked_bindings(check_tool = function(tool) tool, .package = "bambamR")
  local_mocked_bindings(system2 = function(...) structure("bad index", status = 42L),
                        .package = "base")
  expect_error(bambamR:::.align_star("reads.fastq", "index", withr::local_tempdir(),
                                     1L, FALSE, NULL), "STAR.*42")
})

test_that("FASTQ fallback rejects truncation and invalid record structure", {
  fq <- withr::local_tempfile()
  for (lines in list(c("@r", "AC", "+"), c("r", "AC", "+", "II"),
                     c("@r", "AC", "bad", "II"), c("@r", "AC", "+", "I"))) {
    writeLines(lines, fq)
    expect_error(bambamR:::.read_fastq_base(fq, NULL), "FASTQ")
  }
  writeLines(c("@r", "AC", "+", "II", "truncated"), fq)
  expect_equal(nrow(bambamR:::.read_fastq_base(fq, 1L)), 1L)
  expect_error(bb_read_fastq(fq, n = -1), "non-negative")
})

test_that("limma applies a single contrast in a larger design", {
  skip_if_not_installed("limma")
  skip_if_not_installed("edgeR")
  set.seed(123)
  counts <- matrix(rnbinom(1800, mu = 50, size = 5), 200,
                   dimnames = list(paste0("g", 1:200), paste0("s", 1:9)))
  group <- factor(rep(c("a", "b", "c"), each = 3))
  design <- model.matrix(~ group)
  contrast <- c(0, 1, -1)
  result <- bb_limma_voom(counts, design, contrast)
  dge <- edgeR::calcNormFactors(edgeR::DGEList(counts))
  fit <- limma::lmFit(limma::voom(dge, design, plot = FALSE), design)
  fit <- limma::eBayes(limma::contrasts.fit(fit, contrast))
  expected <- limma::topTable(fit, coef = 1, number = nrow(counts), sort.by = "none")
  expect_equal(result$log2fc, expected$logFC)
  expect_equal(result$padj, expected$adj.P.Val)
})

test_that("oncoplot gene selection retains the input cohort denominator", {
  mutations <- data.frame(sample = c("s1", "s2", "s3"),
                          gene = c("A", "B", "B"), mutation_type = "Missense_Mutation")
  plot <- bb_oncoplot(mutations, genes = "A", show_barplot = FALSE)
  expect_equal(plot$scales$get_scales("x")$limits, c("s1", "s2", "s3"))
  expect_equal(unname(plot$scales$get_scales("y")$labels), "A (33%)")
})

test_that("sample and gene separators do not create false multi-hit mutations", {
  mutations <- data.frame(sample = c("a___b", "a"), gene = c("c", "b___c"),
                          mutation_type = "Missense_Mutation")
  expect_identical(bambamR:::.handle_multi_hit(mutations), mutations)
})

test_that("heatmap annotations are displayed and unscaled legends are truthful", {
  counts <- integrity_counts()
  plot <- bb_heatmap(counts, scale = "none")
  expect_equal(plot$scales$get_scales("fill")$name, "Expression")
  skip_if_not_installed("patchwork")
  meta <- data.frame(group = c("a", "b"), row.names = colnames(counts))
  plot <- bb_heatmap(counts, annotation_col = meta)
  expect_s3_class(plot, "patchwork")
  expect_equal(as.character(plot$data$annotation_value), c("a", "b"))
})

test_that("edgeR pipeline preserves adjustment variables in the supplied design", {
  counts <- matrix(rep(1:6, 4), 4,
                   dimnames = list(paste0("g", 1:4), paste0("s", 1:6)))
  meta <- data.frame(condition = rep(c("a", "b"), 3),
                     batch = rep(c("x", "y", "z"), each = 2),
                     row.names = colnames(counts))
  observed <- NULL
  local_mocked_bindings(bb_edger = function(counts, group, design) {
    observed <<- design
    "result"
  }, .package = "bambamR")
  bambamR:::.run_de(counts, meta[6:1, ], "edgeR", ~ batch + condition)
  expect_equal(observed, model.matrix(~ batch + condition, meta))
})
