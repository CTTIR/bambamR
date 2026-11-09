# Launch bambamR Shiny App

Starts the interactive bambamR Shiny application for RNA-seq analysis.
Requires the `shiny` package.

## Usage

``` r
bb_run_app(...)
```

## Arguments

- ...:

  Additional arguments passed to
  [`shiny::runApp()`](https://rdrr.io/pkg/shiny/man/runApp.html).

## Value

This function does not return a value; it launches a Shiny app.

## Examples

``` r
# \donttest{
if (interactive() && requireNamespace("shiny", quietly = TRUE)) {
  bb_run_app()
}
# }
```
