# Print method for bb_qc

Print method for bb_qc

## Usage

``` r
# S3 method for class 'bb_qc'
print(x, ...)
```

## Arguments

- x:

  A bb_qc object.

- ...:

  Additional arguments (ignored).

## Value

`x`, invisibly; called for its side effect of printing a summary.

## Examples

``` r
# \donttest{
# bb_qc() reads FASTQ files with the optional 'ShortRead' package when it
# is installed; loading it can take several seconds.
tmp <- tempfile(fileext = ".fastq")
writeLines(c("@r1", "ACGT", "+", "IIII"), tmp)
qc <- bb_qc(fastq_path = tmp, use_fastqc = FALSE)
print(qc)
#> bambamR QC Summary
#> ==================
#> Files analyzed: 1 
#>   file19af584c280e.fastq: 1 reads
unlink(tmp)
# }
```
