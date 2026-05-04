#' Launch bambamR Shiny App
#'
#' Starts the interactive bambamR Shiny application for RNA-seq analysis.
#' Requires the `shiny` package.
#'
#' @param ... Additional arguments passed to [shiny::runApp()].
#'
#' @return This function does not return a value; it launches a Shiny app.
#'
#' @examples
#' \donttest{
#' if (interactive() && requireNamespace("shiny", quietly = TRUE)) {
#'   bb_run_app()
#' }
#' }
#'
#' @export
bb_run_app <- function(...) {
  check_pkg("shiny", reason = "for the interactive app")

  app_dir <- system.file("app", package = "bambamR")
  if (app_dir == "") {
    stop("Could not find app directory. Try re-installing 'bambamR'.",
         call. = FALSE)
  }
  shiny::runApp(app_dir, ...)
}
