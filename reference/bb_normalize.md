# Normalize Count Matrix

Supports CPM, TPM, TMM (edgeR), and RLE (DESeq2) normalization. CPM and
TPM are always available. TMM requires edgeR, RLE requires DESeq2.

## Usage

``` r
bb_normalize(
  counts,
  method = c("cpm", "tpm", "tmm", "rle"),
  gene_lengths = NULL
)
```

## Arguments

- counts:

  Numeric matrix. Raw count matrix with gene IDs as rownames and sample
  names as colnames.

- method:

  Character. One of `"cpm"`, `"tpm"`, `"tmm"`, `"rle"`.

- gene_lengths:

  Numeric vector. Gene lengths in base pairs. Required for TPM
  normalization. Must be named or in the same order as rownames of
  `counts`.

## Value

A normalized numeric matrix with the same dimensions as `counts`.

## Details

- **CPM** (Counts Per Million): `counts / library_size * 1e6`. Pure base
  R.

- **TPM** (Transcripts Per Million): Normalizes by gene length then
  library size. Pure base R but requires `gene_lengths`.

- **TMM** (Trimmed Mean of M-values): Uses
  [`edgeR::calcNormFactors()`](https://rdrr.io/pkg/edgeR/man/calcNormFactors.html).
  Requires the `edgeR` package.

- **RLE** (Relative Log Expression): Uses
  `DESeq2::estimateSizeFactors()`. Requires the `DESeq2` package.

## Examples

``` r
counts <- matrix(
  rpois(600, lambda = 100),
  nrow = 100, ncol = 6,
  dimnames = list(paste0("gene", 1:100), paste0("sample", 1:6))
)
cpm <- bb_normalize(counts, method = "cpm")
head(cpm)
#>         sample1   sample2   sample3   sample4   sample5   sample6
#> gene1 10118.213  9637.587  8835.501  9311.540 11479.207  8033.323
#> gene2  9617.311 10440.719  9232.602 10797.424  9767.395  8727.561
#> gene3 11019.836 10039.153 11118.833  9509.658  9263.921 10016.860
#> gene4  9216.590  9737.978 10324.630 10599.307 10774.343  9322.622
#> gene5  9917.852 10942.676  9232.602  8717.187  9868.090 11008.628
#> gene6  9216.590 10139.544  9530.428  8519.069  9566.005  6843.201

# TPM requires gene lengths
gene_lengths <- sample(500:5000, 100)
tpm <- bb_normalize(counts, method = "tpm", gene_lengths = gene_lengths)
```
