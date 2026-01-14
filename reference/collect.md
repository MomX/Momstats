# Collect results from statistical analyses

Generic function to extract and add results from statistical analyses to
tibbles. Each method adds appropriate columns based on the analysis
type.

## Usage

``` r
collect(x, ...)
```

## Arguments

- x:

  An object from a statistical analysis (e.g., `stat_pca`, `stat_lda`)

- ...:

  Additional arguments passed to methods

## Value

A tibble with results added

## Details

Different statistical methods add different types of results:

- `stat_pca`: PC scores

- `stat_lda`: Discriminant scores, predictions, posterior probabilities

- `stat_kmeans`: Cluster assignments

## See also

[`collect.stat_pca()`](https://momx.github.io/Momstats/reference/collect.stat_pca.md)
