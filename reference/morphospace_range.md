# Morphospace visualization: data range grid

Create shapes on a grid within the actual data bounding box.

## Usage

``` r
morphospace_range(
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

  Numeric. Extent factor (not used for grid placement). Default 1.0.

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

## Details

Creates a grid within the actual data bounding box. Grid always stays
within data range regardless of extent.

## Examples

``` r
if (FALSE) { # \dontrun{
pca <- boteft %>% stat_pca()
plot(pca, morphospace = "range")
} # }
```
