# Hierarchical clustering for morphometric data

Perform hierarchical clustering on coefficient data.

## Usage

``` r
stat_hclust(
  data,
  formula = NULL,
  method = "ward.D2",
  dist_method = "euclidean",
  center = TRUE,
  scale = NULL,
  k = NULL,
  h = NULL,
  ...
)
```

## Arguments

- data:

  A tibble with coefficient columns

- formula:

  A formula specifying predictors. Can be:

  - Missing: auto-detects single coe column

  - Bare column name: `coe`

  - Formula: `~ coe`, `~ coe + size`, `~ coe1 + coe2`

- method:

  Character. Agglomeration method for hierarchical clustering. One of
  "ward.D2" (default), "single", "complete", "average", "mcquitty",
  "median", "centroid". See
  [`stats::hclust()`](https://rdrr.io/r/stats/hclust.html).

- dist_method:

  Character. Distance metric. Default is "euclidean". Other options:
  "manhattan", "maximum", "canberra", "binary", "minkowski". See
  [`stats::dist()`](https://rdrr.io/r/stats/dist.html).

- center:

  Logical. Should data be centered? Default `TRUE`.

- scale:

  Logical or NULL. Should data be scaled to unit variance? If `NULL`
  (default), automatically determined based on predictor types.

- k:

  Integer. Optional. If provided, cuts tree to k clusters.

- h:

  Numeric. Optional. If provided, cuts tree at height h.

- ...:

  Additional arguments passed to
  [`stats::hclust()`](https://rdrr.io/r/stats/hclust.html)

## Value

An object of class `c("stat_hclust", "momstats")` containing:

- `data`: Original tibble (unchanged)

- `model`: The [`stats::hclust()`](https://rdrr.io/r/stats/hclust.html)
  object

- `dist_matrix`: Distance matrix used for clustering

- `method`: Agglomeration method

- `dist_method`: Distance metric used

- `call`: The function call

- `formula`: Formula used (if any)

- `predictor_cols`: All predictor column names

- `center`: Logical, was centering applied

- `scale`: Logical, was scaling applied

- `k`: Number of clusters (if tree was cut)

- `h`: Height threshold (if tree was cut by height)

- `clusters`: Cluster assignments (if tree was cut)

## Details

`stat_hclust()` provides hierarchical clustering for morphometric data
with proper handling of coefficient columns and optional covariates.

### Agglomeration methods

- `"ward.D2"` (default): Ward's minimum variance method - typically best
  for morphometric data as it minimizes within-cluster variance

- `"complete"`: Maximum distance between clusters

- `"average"`: UPGMA - average distance between clusters

- `"single"`: Minimum distance (tends to chain)

### Distance metrics

- `"euclidean"` (default): Standard L2 distance

- `"manhattan"`: L1 distance, more robust to outliers

- `"maximum"`: Chebyshev distance

### Cutting the tree

The tree can be cut during creation or later:

    # Cut during creation
    hc <- stat_hclust(data, k = 4)

    # Cut later via collect
    hc <- stat_hclust(data)
    data_clustered <- collect(hc, k = 4)

### Getting results

Use [`collect()`](https://momx.github.io/Momstats/reference/collect.md)
to add cluster assignments:

    hc <- boteft %>% stat_hclust()
    boteft_clustered <- collect(hc, k = 4)

Use
[`transduce()`](https://momx.github.io/Momstats/reference/transduce.md)
to reconstruct shapes:

    # Cluster centers (if tree was cut)
    centers <- transduce(hc, tibble(cluster = 1:4))

    # Any node in the tree
    node_shapes <- transduce(hc, tibble(node = c(45, 50, 55)))

## See also

[`stats::hclust()`](https://rdrr.io/r/stats/hclust.html),
[`stats::dist()`](https://rdrr.io/r/stats/dist.html),
[`collect.stat_hclust()`](https://momx.github.io/Momstats/reference/collect.stat_hclust.md),
[`plot.stat_hclust()`](https://momx.github.io/Momstats/reference/plot.stat_hclust.md),
[`transduce.stat_hclust()`](https://momx.github.io/Momstats/reference/transduce.stat_hclust.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic hierarchical clustering
hc1 <- boteft %>% stat_hclust()

# With specific method
hc2 <- boteft %>% stat_hclust(method = "average")

# Cut tree during creation
hc3 <- boteft %>% stat_hclust(k = 4)

# Different distance
hc4 <- boteft %>% stat_hclust(dist_method = "manhattan")

# Add cluster assignments
boteft_clustered <- collect(hc1, k = 4)

# Get cluster center shapes
centers <- transduce(hc3, tibble(cluster = 1:4))

# Get shapes at internal nodes
nodes <- transduce(hc1, tibble(node = c(45, 50)))

# Plot (requires ape package)
plot(hc1)  # Unrooted phylogram (default)
plot(hc1, color = type)  # Color by grouping
plot(hc1, type = "dendrogram")  # Classic dendrogram
} # }
```
