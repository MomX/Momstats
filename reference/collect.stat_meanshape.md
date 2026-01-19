# Collect mean shapes and add to original data

Add mean shapes back to the original data, matched by group.

## Usage

``` r
# S3 method for class 'stat_meanshape'
collect(x, data = NULL, prefix = "mean_", ...)
```

## Arguments

- x:

  A `stat_meanshape` object

- data:

  A tibble. If NULL, uses the original data from the analysis.

- prefix:

  Character. Prefix for mean shape columns. Default "mean\_".

- ...:

  Additional arguments (reserved)

## Value

A tibble with mean shape columns added (and their `_i` inverses)

## Details

For grouped mean shapes, mean shapes are matched back to original rows
by group membership. For global mean shapes, the same mean is added to
all rows.

New columns added:

- `{prefix}{coe_col}`: Mean coefficient column (e.g., `mean_VD`)

- `{prefix}{coe_col}_i`: Mean shape (e.g., `mean_VD_i`)

## Examples

``` r
if (FALSE) { # \dontrun{
ms <- vitis %>% stat_meanshape(status)

# Add mean shapes to data
vitis_with_means <- collect(ms)

# Custom prefix
vitis_with_means <- collect(ms, prefix = "grp_mean_")
} # }
```
