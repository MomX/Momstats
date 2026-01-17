# Morphospace visualization: shapes along axes

Create shapes along specified axes.

## Usage

``` r
morphospace_axes(
  object,
  pcs = c(1, 2),
  n = 1,
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

- n:

  Integer. Number of shapes per axis direction. Default 1.

- extent:

  Numeric. Extent factor (\< 1 = zoom in, \> 1 = extrapolate). Default
  1.0.

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

Creates shapes along each axis. For n=1, creates 4 shapes at extremes.
For n\>1, creates more intermediate shapes along each axis.

## Examples

``` r
if (FALSE) { # \dontrun{
pca <- boteft %>% stat_pca()
plot(pca, morphospace = "axes")

# More shapes
plot(pca, morphospace = FALSE)
morphospace_axes(pca, n = 3)
} # }
```
