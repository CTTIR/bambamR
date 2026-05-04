# Script to generate example_counts.rds
# Run once to create the bundled example dataset
set.seed(42)

n_genes <- 200
n_samples <- 10

# Simulate count data with some DE genes
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

rownames(counts) <- paste0("GENE", seq_len(n_genes))
colnames(counts) <- paste0("S", seq_len(n_samples))

metadata <- data.frame(
  condition = factor(rep(c("control", "treatment"), each = 5)),
  batch = factor(rep(c("A", "B"), times = 5)),
  row.names = colnames(counts)
)

example_data <- list(counts = counts, metadata = metadata)
saveRDS(example_data, "inst/extdata/example_counts.rds")
