# Test helper: build a tiny sorted, indexed BAM file from an inline SAM record
# set. Returns the BAM path, or NULL if Rsamtools is unavailable so callers can
# skip. The BAM is created inside `dir` (default a fresh temp dir).

make_test_bam <- function(dir = tempfile("bb_bam"),
                          n_mapped = 2L, n_unmapped = 1L) {
  if (!requireNamespace("Rsamtools", quietly = TRUE)) {
    return(NULL)
  }
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)

  seq50 <- paste(rep("A", 50), collapse = "")
  qual50 <- paste(rep("I", 50), collapse = "")

  records <- character(0)
  for (i in seq_len(n_mapped)) {
    pos <- 100L + (i - 1L) * 200L
    records <- c(records, paste(
      paste0("r", i), "0", "chr1", pos, "60", "50M", "*", "0", "0",
      seq50, qual50, sep = "\t"
    ))
  }
  for (j in seq_len(n_unmapped)) {
    records <- c(records, paste(
      paste0("u", j), "4", "*", "0", "0", "*", "*", "0", "0",
      seq50, qual50, sep = "\t"
    ))
  }

  sam <- file.path(dir, "reads.sam")
  writeLines(c("@HD\tVN:1.6\tSO:coordinate", "@SQ\tSN:chr1\tLN:2000",
               records), sam)

  tryCatch(
    Rsamtools::asBam(sam, destination = file.path(dir, "reads"),
                     overwrite = TRUE),
    error = function(e) NULL
  )
}
