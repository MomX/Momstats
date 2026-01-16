# Plot k-means clustering results

Visualize k-means clustering with various plot types.

## Usage

``` r
# S3 method for class 'stat_kmeans'
plot(x, type = c("clusters", "withinss"), color = NULL, labels = NULL, ...)
```

## Arguments

- x:

  A `stat_kmeans` object

- type:

  Character. Type of plot:

  - `"clusters"`: Cluster visualization (requires PCA projection)

  - `"withinss"`: Within-cluster sum of squares by cluster

- color:

  Column name (bare or quoted) for coloring points (for cluster plot).
  If NULL, colors by cluster assignment.

- labels:

  Column name (bare or quoted) for text labels.

- ...:

  Additional arguments (reserved)

## Value

NULL (invisibly). Draws plot as side effect.

## Examples

``` r
if (FALSE) { # \dontrun{
km <- boteft %>% stat_kmeans(k = 3)

# Cluster plot (colored by cluster)
plot(km)

# Color by original grouping
plot(km, color = type)

# Within-cluster SS plot
plot(km, type = "withinss")
} # }
```
