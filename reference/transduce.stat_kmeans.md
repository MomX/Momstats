# Transduce k-means cluster centers to shapes

Reconstruct shapes at cluster centers or at positions relative to
centers.

## Usage

``` r
# S3 method for class 'stat_kmeans'
transduce(object, positions)
```

## Arguments

- object:

  A `stat_kmeans` object

- positions:

  A tibble with cluster specifications. Must contain a `cluster` column
  with cluster numbers (1 to k).

## Value

A tibble with reconstructed coefficients and shapes

## Details

For k-means, transduce simply retrieves the cluster centers from the
model and reconstructs the corresponding shapes. Since cluster centers
are already in the original predictor space, no inverse transformation
is needed.

## Examples

``` r
if (FALSE) { # \dontrun{
km <- boteft %>% stat_kmeans(k = 3)

# Get shapes at all cluster centers
centers <- transduce(km, tibble(cluster = 1:3))

# Get shape at specific cluster
center1 <- transduce(km, tibble(cluster = 1))
} # }
```
