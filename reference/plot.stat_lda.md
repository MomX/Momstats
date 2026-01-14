# Plot LDA results

Visualize LDA results with discrimination plots, loading plots, or
histograms (for binary classification).

## Usage

``` r
# S3 method for class 'stat_lda'
plot(
  x,
  type = c("scores", "loadings"),
  lds = c(1, 2),
  color = NULL,
  labels = NULL,
  chull = TRUE,
  cex = 0.7,
  legend = TRUE,
  ...
)
```

## Arguments

- x:

  A `stat_lda` object

- type:

  Character. Type of plot:

  - `"scores"`: LD score plot (default) - scatter for 2+ LDs, boxplot
    for 1 LD

  - `"loadings"`: Loading plot for first two LDs

- lds:

  Integer vector of length 2. Which LDs to plot for score/loading plots.
  Default is `c(1, 2)`. For binary classification (1 LD), only first
  element used.

- color:

  Column name (bare or quoted) for coloring points. Default uses the
  response variable.

- labels:

  Column name (bare or quoted) for text labels in score plot. For
  loadings, logical: should variable names be shown? Default `TRUE`.

- chull:

  Logical. Draw convex hulls around groups? Only works when `color` is a
  factor. Default `TRUE`.

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
lda <- boteft %>% stat_lda(type)

# Score plot (default)
plot(lda)
plot(lda, color = type)
plot(lda, labels = pred)
plot(lda, color = type, chull = TRUE)

# Different LDs
plot(lda, lds = c(2, 3))

# Loading plot
plot(lda, type = "loadings")
plot(lda, type = "loadings", labels = FALSE)
} # }
```
