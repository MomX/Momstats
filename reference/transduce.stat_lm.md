# Transduce linear model positions to shapes

Predict coefficients and shapes at specified predictor values.

## Usage

``` r
# S3 method for class 'stat_lm'
transduce(object, positions)
```

## Arguments

- object:

  A `stat_lm` object

- positions:

  A tibble with predictor values. Column names must match the model's
  predictors.

## Value

A tibble with predicted values

## Examples

``` r
if (FALSE) { # \dontrun{
lm_fit <- boteft %>% stat_lm(coe ~ length)

# Predict at new lengths
new_shapes <- transduce(lm_fit, tibble::tibble(length = seq(50, 150, by = 10)))

# Multiple predictors
lm_fit2 <- boteft %>% stat_lm(coe ~ length + width)
new_shapes <- transduce(lm_fit2,
                        expand_grid(length = c(50, 100, 150),
                                    width = c(20, 30)))
} # }
```
