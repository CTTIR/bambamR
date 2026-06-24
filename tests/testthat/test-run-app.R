# Tests for R/run_app.R. We never actually launch the Shiny server; we verify
# the dependency guard and that the bundled app directory resolves.

test_that("bb_run_app errors when shiny is unavailable", {
  skip_if(requireNamespace("shiny", quietly = TRUE),
          "shiny is installed; cannot test the missing-package guard")
  expect_error(bb_run_app(), "shiny")
})

test_that("the bundled Shiny app directory is installed", {
  app_dir <- system.file("app", package = "bambamR")
  expect_true(nzchar(app_dir))
  expect_true(file.exists(file.path(app_dir, "app.R")))
})

test_that("bb_run_app dispatches to shiny::runApp when shiny is present", {
  skip_if_not_installed("shiny")
  called <- new.env()
  called$dir <- NULL
  local_mocked_bindings(
    runApp = function(appDir, ...) {
      called$dir <- appDir
      invisible(NULL)
    },
    .package = "shiny"
  )
  bb_run_app()
  expect_true(nzchar(called$dir))
  expect_true(grepl("app$", called$dir))
})
