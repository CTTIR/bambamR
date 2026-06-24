#!/usr/bin/env Rscript
# fingerprint.R --------------------------------------------------------------
# Deterministic SHA-256 fingerprint of the package's R/ source tree.
#
# Purpose: the coverage work in this package is *test-only* and must preserve
# behaviour. This script produces a content hash of every file under R/ so we
# can prove the implementation is byte-for-byte unchanged before and after the
# test additions.
#
# Usage:
#   Rscript _coverage/fingerprint.R            # write _coverage/fingerprint.txt
#   Rscript _coverage/fingerprint.R --compare  # compare against the saved file
#
# Determinism: files are sorted by path, read as raw bytes, and digested with
# digest::digest(algo = "sha256"). No timestamps or environment state enter the
# hash, so repeated runs on identical sources yield identical output.

suppressWarnings(suppressMessages({
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("The 'digest' package (dev-only) is required to run fingerprint.R.")
  }
}))

# Resolve the R/ directory relative to this script regardless of CWD.
args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
script_dir <- if (length(file_arg) == 1L) {
  dirname(normalizePath(file_arg))
} else {
  getwd()
}
r_dir <- normalizePath(file.path(script_dir, "..", "R"), mustWork = TRUE)
out_file <- file.path(script_dir, "fingerprint.txt")

fingerprint_dir <- function(dir) {
  files <- list.files(dir, full.names = TRUE, recursive = TRUE)
  files <- sort(files[file.info(files)$isdir %in% FALSE])
  rel <- sub(paste0("^", normalizePath(dir), "[\\/]"), "",
             normalizePath(files))
  rel <- gsub("\\\\", "/", rel)
  hashes <- vapply(files, function(f) {
    digest::digest(readBin(f, what = "raw", n = file.info(f)$size),
                   algo = "sha256", serialize = FALSE)
  }, character(1))
  out <- paste0(hashes, "  ", rel)
  # Append an overall hash of the per-file lines for a single comparison token.
  overall <- digest::digest(paste(out, collapse = "\n"),
                            algo = "sha256", serialize = FALSE)
  c(out, paste0(overall, "  <OVERALL>"))
}

trailing <- commandArgs(trailingOnly = TRUE)
compare <- "--compare" %in% trailing

current <- fingerprint_dir(r_dir)

if (compare) {
  if (!file.exists(out_file)) {
    stop("No saved fingerprint to compare against: ", out_file)
  }
  saved <- readLines(out_file)
  if (identical(saved, current)) {
    cat("FINGERPRINT: IDENTICAL\n")
  } else {
    cat("FINGERPRINT: CHANGED\n")
    added <- setdiff(current, saved)
    removed <- setdiff(saved, current)
    if (length(removed)) cat("-- only in saved --\n", paste(removed, collapse = "\n"), "\n")
    if (length(added)) cat("-- only in current --\n", paste(added, collapse = "\n"), "\n")
    quit(status = 1L)
  }
} else {
  writeLines(current, out_file)
  cat("Wrote fingerprint for", length(current) - 1L, "R/ files to", out_file, "\n")
  cat("OVERALL:", sub("  <OVERALL>$", "", current[length(current)]), "\n")
}
