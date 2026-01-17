# Morphospace visualization: uniform grid

Create shapes on a uniform grid across the full PC plane.

## Usage

``` r
morphospace_grid(
  object,
  pcs = c(1, 2),
  nrow = 5,
  ncol = 5,
  extent = 1,
  transduce = TRUE,
  template = TRUE,
  draw = TRUE,
  ...
)
```

## Arguments

- object:

  A stat\_\* object (stat_pca, stat_lda, etc.)

- pcs:

  Integer vector of length 2. Which axes to use. Default c(1, 2).

- nrow:

  Integer. Number of rows in grid. Default 5.

- ncol:

  Integer. Number of columns in grid. Default 5.

- extent:

  Numeric. Extent factor. Default 1.0.

- transduce:

  Logical. Should shapes be transduced? Default TRUE.

- template:

  Logical or numeric. Templating option. Default TRUE.

- draw:

  Logical. Should shapes be drawn? Default TRUE.

- ...:

  Additional graphical parameters (lwd, col, etc.)

## Value

A tibble with positions, coefficients, and shapes

## Examples

``` r
if (FALSE) { # \dontrun{
pca <- boteft %>% stat_pca()
plot(pca, morphospace = "grid", nrow = 3, ncol = 3)
} # }
```
