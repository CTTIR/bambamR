# Compatibility entry point for the registered v5 package artwork.
# Usage: Rscript data-raw/create_logo.R

args <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", args[grepl("^--file=", args)])
root <- if (length(script)) {
  normalizePath(file.path(dirname(script[1]), ".."), mustWork = TRUE)
} else {
  normalizePath(".", mustWork = TRUE)
}
renderer <- file.path(root, "tools", "branding", "render.R")
if (!file.exists(renderer)) stop("Run this entry point from the package source tree.")
status <- system2(file.path(R.home("bin"), "Rscript"),
                  c(shQuote(renderer), "--apply"))
if (!identical(status, 0L)) stop("Package artwork regeneration failed.")
