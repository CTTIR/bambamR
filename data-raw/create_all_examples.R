# Generate all example datasets for bambamR
# Run this script once to create bundled example data
# All datasets must be < 1MB total

set.seed(42)

# ============================================================
# 1. RNA-seq count matrix + metadata (existing, enhanced)
# ============================================================
n_genes <- 200
n_samples <- 10

# Realistic gene names (mix of real gene symbols)
gene_names <- c(
  paste0("TP53"), "KRAS", "EGFR", "BRAF", "PIK3CA",
  "PTEN", "APC", "CDKN2A", "RB1", "MYC",
  "BRCA1", "BRCA2", "ATM", "NRAS", "IDH1",
  "SMAD4", "VHL", "NF1", "NOTCH1", "FBXW7",
  paste0("GENE", 21:200)
)

base_means <- rpois(n_genes, lambda = 500)
counts <- matrix(0L, nrow = n_genes, ncol = n_samples)
for (j in seq_len(n_samples)) {
  counts[, j] <- rpois(n_genes, lambda = base_means)
}

# Add differential expression for first 20 genes in treatment group
for (j in 6:10) {
  counts[1:10, j] <- rpois(10, lambda = base_means[1:10] * 3)
  counts[11:20, j] <- rpois(10, lambda = base_means[11:20] * 0.3)
}

rownames(counts) <- gene_names
colnames(counts) <- paste0("Sample_", sprintf("%02d", seq_len(n_samples)))

metadata <- data.frame(
  condition = factor(rep(c("control", "treatment"), each = 5)),
  batch = factor(rep(c("A", "B"), times = 5)),
  sex = factor(c("M", "F", "M", "F", "M", "F", "M", "F", "M", "F")),
  age = c(45, 52, 38, 61, 55, 47, 59, 42, 50, 44),
  row.names = colnames(counts)
)

example_counts <- list(counts = counts, metadata = metadata)
saveRDS(example_counts, "inst/extdata/example_counts.rds")
cat("example_counts.rds:", file.size("inst/extdata/example_counts.rds"), "bytes\n")


# ============================================================
# 2. Mutation data for oncoplot
# ============================================================
n_mut_samples <- 50
n_mutations <- 300

mut_types <- c("Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del",
               "Frame_Shift_Ins", "Splice_Site", "In_Frame_Del",
               "In_Frame_Ins", "Translation_Start_Site")

onco_genes <- c("TP53", "KRAS", "PIK3CA", "PTEN", "APC",
                "BRAF", "EGFR", "NRAS", "CDKN2A", "RB1",
                "SMAD4", "ATM", "BRCA1", "BRCA2", "IDH1",
                "VHL", "NF1", "NOTCH1", "FBXW7", "ARID1A")

# Weight genes by realistic mutation frequency
gene_weights <- c(0.18, 0.12, 0.10, 0.08, 0.08,
                  0.06, 0.06, 0.05, 0.05, 0.04,
                  0.03, 0.03, 0.02, 0.02, 0.02,
                  0.01, 0.01, 0.01, 0.01, 0.02)

mut_data <- data.frame(
  sample = sample(paste0("TCGA-", sprintf("%03d", seq_len(n_mut_samples))),
                  n_mutations, replace = TRUE),
  gene = sample(onco_genes, n_mutations, replace = TRUE, prob = gene_weights),
  mutation_type = sample(mut_types, n_mutations, replace = TRUE,
                         prob = c(0.45, 0.15, 0.10, 0.05, 0.10,
                                  0.05, 0.03, 0.07)),
  stringsAsFactors = FALSE
)

# Clinical annotation for the samples
all_samples <- unique(mut_data$sample)
clinical_data <- data.frame(
  Stage = sample(c("I", "II", "III", "IV"), length(all_samples),
                 replace = TRUE, prob = c(0.15, 0.30, 0.35, 0.20)),
  Gender = sample(c("Male", "Female"), length(all_samples), replace = TRUE),
  Smoking = sample(c("Never", "Former", "Current"), length(all_samples),
                   replace = TRUE, prob = c(0.4, 0.35, 0.25)),
  row.names = all_samples
)

example_mutations <- list(mutations = mut_data, clinical = clinical_data)
saveRDS(example_mutations, "inst/extdata/example_mutations.rds")
cat("example_mutations.rds:", file.size("inst/extdata/example_mutations.rds"), "bytes\n")


# ============================================================
# 3. Pre-computed DE results (for users who just want to plot)
# ============================================================
n_de <- 500
de_results <- data.frame(
  gene = c(gene_names, paste0("GENE", 201:500)),
  log2fc = c(rnorm(20, mean = c(rep(2, 10), rep(-2, 10)), sd = 0.5),
             rnorm(480, mean = 0, sd = 0.5)),
  pvalue = c(runif(20, 1e-10, 1e-3), runif(480, 0.01, 1)),
  padj = NA_real_,
  basemean = 10^runif(500, 1, 4),
  stringsAsFactors = FALSE
)
de_results$padj <- p.adjust(de_results$pvalue, method = "BH")

example_de <- de_results
saveRDS(example_de, "inst/extdata/example_de_results.rds")
cat("example_de_results.rds:", file.size("inst/extdata/example_de_results.rds"), "bytes\n")


# ============================================================
# 4. Small FASTQ file (for testing import)
# ============================================================
n_reads <- 100
bases <- c("A", "C", "G", "T")
read_len <- 75

fastq_lines <- character(n_reads * 4)
for (i in seq_len(n_reads)) {
  idx <- (i - 1) * 4
  seq_str <- paste0(sample(bases, read_len, replace = TRUE), collapse = "")
  qual_str <- paste0(sample(c("I", "H", "G", "F", "E", "D", "C", "B", "A",
                               "?", ">", "="), read_len, replace = TRUE,
                            prob = c(0.20, 0.18, 0.15, 0.12, 0.10, 0.08,
                                     0.06, 0.04, 0.03, 0.02, 0.01, 0.01)),
                     collapse = "")
  fastq_lines[idx + 1] <- paste0("@read_", sprintf("%04d", i),
                                  " length=", read_len)
  fastq_lines[idx + 2] <- seq_str
  fastq_lines[idx + 3] <- "+"
  fastq_lines[idx + 4] <- qual_str
}
writeLines(fastq_lines, "inst/extdata/example_reads.fastq")
cat("example_reads.fastq:", file.size("inst/extdata/example_reads.fastq"), "bytes\n")


# ============================================================
# 5. Gene lengths (for TPM normalization example)
# ============================================================
gene_lengths <- data.frame(
  gene = gene_names,
  length = sample(300:15000, n_genes, replace = TRUE),
  stringsAsFactors = FALSE
)
saveRDS(gene_lengths, "inst/extdata/example_gene_lengths.rds")
cat("example_gene_lengths.rds:", file.size("inst/extdata/example_gene_lengths.rds"), "bytes\n")

cat("\nTotal extdata size:",
    sum(file.size(list.files("inst/extdata", full.names = TRUE,
                             pattern = "\\.(rds|fastq)$"))),
    "bytes\n")
