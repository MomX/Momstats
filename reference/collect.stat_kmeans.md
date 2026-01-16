# Collect cluster assignments from k-means

Extract cluster assignments from k-means results and add them to a
tibble.

## Usage

``` r
# S3 method for class 'stat_kmeans'
collect(x, data = NULL, name = "cluster", ...)
```

## Arguments

- x:

  A `stat_kmeans` object

- data:

  A tibble. If NULL, uses the original data from the clustering.

- name:

  Character. Name for the cluster assignment column. Default is
  "cluster".

- ...:

  Additional arguments (reserved)

## Value

A tibble with cluster assignments added

## Examples

``` r
if (FALSE) { # \dontrun{
km <- boteft %>% stat_kmeans(k = 3)

# Add cluster assignments
collect(km)

# Custom column name
collect(km, name = "k3_cluster")
} # }
```
