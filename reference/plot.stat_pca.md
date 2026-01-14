# Plot PCA results

Visualize PCA results with score plots, scree plots, or loading plots.

## Usage

``` r
# S3 method for class 'stat_pca'
plot(
  x,
  type = c("scores", "scree", "loadings"),
  pcs = c(1, 2),
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

  Character. Type of plot:

  - `"scores"`: PC score plot (default)

  - `"scree"`: Variance explained by each PC

  - `"loadings"`: Loading plot for first two PCs

- pcs:

  Integer vector of length 2. Which PCs to plot for score/loading plots.
  Default is `c(1, 2)`.

- color:

  Column name (bare or quoted) for coloring points in score plot.

- labels:

  Column name (bare or quoted) for text labels in score plot. For
  loadings, logical: should variable names be shown? Default `TRUE`.

- chull:

  Logical. Draw convex hulls around groups? Only works when `color` is a
  factor. Default `TRUE`.

- n_pcs:

  Integer. For scree plot, how many PCs to show? Default is all.

- cex:

  Numeric. Character expansion for labels. Default is 0.7.

- legend:

  Logical. Show legend for colors? Default `TRUE`.

- ...:

  Additional arguments (reserved for future use)

## Value

NULL (invisibly). Draws plot as side effect.

## Examples

``` r
if (FALSE) { # \dontrun{
pca <- boteft %>% stat_pca()

# Score plot (default)
plot(pca)
plot(pca, color = type)
plot(pca, labels = type, color = type)
plot(pca, color = type, chull = TRUE)

# Different PCs
plot(pca, pcs = c(2, 3))

# Scree plot
plot(pca, type = "scree")
plot(pca, type = "scree", n_pcs = 10)  # First 10 only

# Loading plot
plot(pca, type = "loadings")
plot(pca, type = "loadings", labels = FALSE)
} # }
```
