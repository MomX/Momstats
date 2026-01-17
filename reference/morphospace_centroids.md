# Morphospace visualization: group centroids

Create shapes at group centroid positions.

## Usage

``` r
morphospace_centroids(
  object,
  group = NULL,
  group_vals = NULL,
  transduce = TRUE,
  template = TRUE,
  draw = TRUE,
  ...
)
```

## Arguments

- object:

  A stat\_\* object (stat_pca, stat_lda, etc.)

- group:

  Column name for grouping (bare or quoted). For PCA, required. For LDA,
  optional (uses discriminant groups if NULL).

- group_vals:

  Pre-evaluated group values (used internally by plot).

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
plot(pca, color = type, morphospace = "centroids")

# Or call directly
plot(pca, morphospace = FALSE)
morphospace_centroids(pca, group = type)
} # }
```
