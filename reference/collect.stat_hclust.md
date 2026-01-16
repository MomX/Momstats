# Collect cluster assignments from hierarchical clustering

Extract cluster assignments from hierarchical clustering and add them to
a tibble.

## Usage

``` r
# S3 method for class 'stat_hclust'
collect(x, data = NULL, k = NULL, h = NULL, name = "cluster", ...)
```

## Arguments

- x:

  A `stat_hclust` object

- data:

  A tibble. If NULL, uses the original data.

- k:

  Integer. Number of clusters to cut tree into.

- h:

  Numeric. Height at which to cut tree.

- name:

  Character. Name for cluster assignment column. Default is "cluster".

- ...:

  Additional arguments (reserved)

## Value

A tibble with cluster assignments added

## Details

If the tree was already cut during
[`stat_hclust()`](https://momx.github.io/Momstats/reference/stat_hclust.md),
those clusters are used unless `k` or `h` is specified here.

## Examples

``` r
if (FALSE) { # \dontrun{
hc <- boteft %>% stat_hclust()

# Cut and collect
collect(hc, k = 4)

# Different k
collect(hc, k = 5)

# Cut by height
collect(hc, h = 10)

# Custom column name
collect(hc, k = 3, name = "hc_cluster")
} # }
```
