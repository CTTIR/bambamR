# Refresh shared branding without rebuilding documentation or examples.
sync_site_branding <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  script <- sub("^--file=", "", args[grepl("^--file=", args)])
  if (length(script) != 1L) stop("Run this file with Rscript.")
  branding <- dirname(normalizePath(script, winslash = "/", mustWork = TRUE))
  repository <- normalizePath(file.path(branding, "..", ".."),
                              winslash = "/", mustWork = TRUE)
  params <- commandArgs(trailingOnly = TRUE)
  if (length(params) != 1L || is.na(params) || !nzchar(params) ||
      !dir.exists(params)) {
    stop("Supply an existing site output directory.")
  }
  output <- normalizePath(params, winslash = "/", mustWork = TRUE)
  figures <- file.path(repository, "man", "figures")
  favicons <- file.path(repository, "pkgdown", "favicon")
  logos <- list.files(figures, "^logo.*[.](svg|png|pdf)$", full.names = TRUE)
  icons <- list.files(favicons, "[.](svg|png|ico|webmanifest)$", full.names = TRUE)
  if (!file.exists(file.path(figures, "logo.svg")) || !length(icons)) {
    stop("The source logo or favicon bundle is missing.")
  }
  roots <- c(output, if (dir.exists(file.path(output, "dev"))) file.path(output, "dev"))
  sources <- destinations <- character()
  for (root in roots) {
    sources <- c(sources, logos, file.path(figures, "logo.svg"), icons)
    destinations <- c(destinations,
      file.path(root, "reference", "figures", basename(logos)),
      file.path(root, "logo.svg"), file.path(root, basename(icons)))
    # Retain and refresh the older raster alias only where it is already used.
    if (file.exists(file.path(root, "logo.png"))) {
      sources <- c(sources, file.path(figures, "logo.png"))
      destinations <- c(destinations, file.path(root, "logo.png"))
    }
  }
  if (!all(file.exists(sources)) || any(file.info(sources)$isdir)) {
    stop("A branding source file is missing or is a directory.")
  }
  for (target in destinations) {
    link <- Sys.readlink(target)
    if (dir.exists(target) || (!is.na(link) && nzchar(link))) {
      stop("A branding destination must be a regular file: ", target)
    }
    ancestor <- dirname(target)
    while (!dir.exists(ancestor)) ancestor <- dirname(ancestor)
    resolved <- normalizePath(ancestor, winslash = "/", mustWork = TRUE)
    if (!identical(resolved, output) && !startsWith(resolved, paste0(output, "/"))) {
      stop("A branding destination escapes the site directory: ", target)
    }
  }
  for (i in seq_along(sources)) {
    dir.create(dirname(destinations[i]), recursive = TRUE, showWarnings = FALSE)
    if (!file.copy(sources[i], destinations[i], overwrite = TRUE)) {
      stop("Could not copy branding asset: ", destinations[i])
    }
  }
  message("Synchronized ", length(destinations), " branding assets across ",
          length(roots), " site root(s).")
}

sync_site_branding()
