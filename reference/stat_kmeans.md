# K-means clustering for morphometric data

Perform k-means clustering on coefficient data.

## Usage

``` r
stat_kmeans(data, formula = NULL, k = 3, nstart = 25, ...)
```

## Arguments

- data:

  A tibble with coefficient columns

- formula:

  A formula specifying predictors. Can be:

  - Missing: auto-detects single coe column

  - Bare column name: `coe`

  - Formula: `~ coe`, `~ coe + size`, `~ coe1 + coe2`

  - Use `coe` in formula to auto-detect coefficient columns

- k:

  Integer. Number of clusters. Default is 3.

- nstart:

  Integer. Number of random starts for k-means. Default is 25.

- ...:

  Additional arguments passed to
  [`stats::kmeans()`](https://rdrr.io/r/stats/kmeans.html)

## Value

An object of class `c("stat_kmeans", "momstats")` containing:

- `data`: Original tibble (unchanged)

- `model`: The [`stats::kmeans()`](https://rdrr.io/r/stats/kmeans.html)
  object

- `method`: "kmeans"

- `call`: The function call

- `formula`: Formula used (if any)

- `predictor_cols`: All predictor column names

- `k`: Number of clusters

- `centers`: Cluster centers in predictor space (matrix)

- `cluster_sizes`: Size of each cluster

- `withinss`: Within-cluster sum of squares

- `tot_withinss`: Total within-cluster sum of squares

- `betweenss`: Between-cluster sum of squares

- `variance_explained`: Proportion of variance explained by clustering

## Details

`stat_kmeans()` provides a unified interface for k-means clustering on
morphometric coefficient data.

### Formula syntax

The formula specifies which predictors to use:

- `~ coe`: Use auto-detected coefficient column(s)

- `~ coe1 + coe2`: Use specific coefficient columns

- `~ coe + size`: Coefficient column plus a covariate

- Bare name or missing: auto-detect single coe column

### Getting results

Use [`collect()`](https://momx.github.io/Momstats/reference/collect.md)
to add cluster assignments to your data:

    km <- boteft %>% stat_kmeans(k = 3)
    boteft_clustered <- collect(km)  # Adds 'cluster' column

Use
[`transduce()`](https://momx.github.io/Momstats/reference/transduce.md)
to get shapes at cluster centers:

    center_shapes <- transduce(km, tibble(cluster = 1:3))

## See also

[`stats::kmeans()`](https://rdrr.io/r/stats/kmeans.html),
[`collect.stat_kmeans()`](https://momx.github.io/Momstats/reference/collect.stat_kmeans.md),
[`plot.stat_kmeans()`](https://momx.github.io/Momstats/reference/plot.stat_kmeans.md),
[`transduce.stat_kmeans()`](https://momx.github.io/Momstats/reference/transduce.stat_kmeans.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic k-means with 3 clusters
km1 <- boteft %>% stat_kmeans(k = 3)

# More clusters
km2 <- boteft %>% stat_kmeans(k = 5)

# With covariate
km3 <- boteft %>% stat_kmeans(~ coe + length, k = 4)

# Add cluster assignments
boteft_clustered <- collect(km1)

# Get cluster center shapes
centers <- transduce(km1, tibble(cluster = 1:3))

# Plot results
plot(km1)  # Cluster visualization
plot(km1, color = type)  # Color by original grouping
} # }
```
