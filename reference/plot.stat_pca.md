# Plot PCA results

Plot PCA results

## Usage

``` r
# S3 method for class 'stat_pca'
plot(
  x,
  type = c("scores", "scree", "loadings"),
  pcs = c(1, 2),
  extent = 1,
  morphospace = "grid",
  color = NULL,
  labels = NULL,
  chull = TRUE,
  n_pcs = NULL,
  cex = 0.7,
  legend = TRUE,
  ...
)
```

## Arguments

- x:

  A `stat_pca` object

- type:

  Character. Type of plot

- pcs:

  Integer vector of length 2. Which PCs to plot

- extent:

  Numeric. Plot extent factor (\< 1 = zoom in, \> 1 = zoom out)

- morphospace:

  Character or FALSE. Morphospace type: "grid", "range", "axes",
  "centroids", FALSE

- color:

  Column name for coloring points

- labels:

  Column name for text labels

- chull:

  Logical. Draw convex hulls?

- n_pcs:

  Integer. For scree plot

- cex:

  Numeric. Character expansion

- legend:

  Logical. Show legend?

- ...:

  Passed to morphospace functions (n, nrow, ncol, template, lwd, col)
