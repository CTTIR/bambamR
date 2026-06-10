# Export Results to RDS

Saves a bambamR result object or any R object as an RDS file.

## Usage

``` r
bb_export_rds(result, path)
```

## Arguments

- result:

  An R object to save (typically a `bb_result` or DE result).

- path:

  Character. File path for the output `.rds` file.

## Value

The `path` invisibly.

## Examples

``` r
tmp <- tempfile(fileext = ".rds")
counts <- matrix(rpois(100, 50), nrow = 10,
  dimnames = list(paste0("g", 1:10), paste0("s", 1:10)))
bb_export_rds(counts, tmp)
#> Saved to: /tmp/Rtmp11ZIBD/file1b962f0ec63.rds
readRDS(tmp)
#>     s1 s2 s3 s4 s5 s6 s7 s8 s9 s10
#> g1  46 61 50 50 54 54 45 60 57  55
#> g2  61 39 41 52 56 49 41 58 60  63
#> g3  48 47 46 53 52 41 49 51 49  42
#> g4  52 52 52 50 50 45 53 47 43  48
#> g5  50 56 53 38 55 38 62 47 54  62
#> g6  59 45 43 48 48 56 47 63 67  49
#> g7  56 51 47 42 44 40 59 62 43  50
#> g8  59 58 43 57 68 41 54 40 44  52
#> g9  52 52 45 53 62 42 46 53 50  52
#> g10 41 49 46 53 41 48 53 49 52  51
```
