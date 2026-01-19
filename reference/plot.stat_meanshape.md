# Plot mean shapes

Visualize mean shapes, optionally with contributing data.

## Usage

``` r
# S3 method for class 'stat_meanshape'
plot(
  x,
  all = FALSE,
  facet = TRUE,
  col = "red",
  lwd = 2,
  all_col = "#000000E6",
  all_lwd = 0.1,
  labels = TRUE,
  label_cex = 0.7,
  ...
)
```

## Arguments

- x:

  A `stat_meanshape` object

- all:

  Logical. Should all contributing shapes be plotted behind the mean?
  Default FALSE.

- facet:

  Logical. For grouped data, should groups be offset vertically? Default
  TRUE. If FALSE, groups are overlaid with colors + legend.

- col:

  Character. Color(s) for mean shapes. Default "red". For facet=FALSE,
  vector of colors per group can be provided.

- lwd:

  Numeric. Line width for mean shapes (or cex for ldk). Default 2.

- all_col:

  Character. Color for contributing shapes when facet=TRUE. Default
  "gray". When facet=FALSE, ignored (uses parent colors with
  transparency).

- all_lwd:

  Numeric. Line width for contributing shapes (or cex for ldk). Default
  0.5.

- labels:

  Logical. Show group labels/legend? Default TRUE.

- label_cex:

  Numeric. Label size. Default 0.7.

- ...:

  Additional graphical parameters (reserved)
