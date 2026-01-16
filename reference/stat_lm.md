# Linear Model for morphometric data

Fit linear models with coefficient data as response or predictor.

## Usage

``` r
stat_lm(data, formula, ...)
```

## Arguments

- data:

  A tibble with coefficient and/or scalar columns

- formula:

  A formula specifying the model. Can be:

  - Shape regression: `coe ~ length` (shape depends on predictors)

  - Response regression: `length ~ coe` (response depends on shape)

  - Multivariate response: `cbind(length, width) ~ coe`

  - Multiple predictors: `coe ~ length + temperature`

  - Interactions: `coe ~ length * temperature`

- ...:

  Additional arguments passed to
  [`stats::lm()`](https://rdrr.io/r/stats/lm.html)

## Value

An object of class `c("stat_lm", "momstats")` containing:

- `data`: Original tibble (unchanged)

- `model`: The [`stats::lm()`](https://rdrr.io/r/stats/lm.html) object

- `method`: "lm"

- `direction`: "shape_regression" or "response_regression"

- `call`: The function call

- `formula`: Formula used

- `response_col`: Response column name(s)

- `predictor_cols`: Predictor column name(s)

- `r_squared`: R^2 value(s)

- `adj_r_squared`: Adjusted R^2 value(s)

## Details

`stat_lm()` provides a unified interface for linear modeling with
morphometric data.

### Formula directions

The function automatically detects model direction from the formula:

**Shape regression** (`coe ~ predictors`):

- Fits multivariate regression with all harmonics as response

- Proper treatment of shape as correlated multivariate outcome

- Tests whether predictors affect shape

- Example: `boteft %>% stat_lm(coe ~ length)`

**Response regression** (`response ~ coe`):

- Fits regression with unfolded harmonics as predictors

- Tests whether shape predicts a response variable

- Example: `boteft %>% stat_lm(length ~ coe)`

**Multivariate response** (`cbind(y1, y2) ~ coe`):

- Fits multivariate regression with multiple responses

- Example: `boteft %>% stat_lm(cbind(length, width) ~ coe)`

### Interactions and multiple predictors

Standard R formula syntax is supported:

- `coe ~ length + temperature`: Additive effects

- `coe ~ length * temperature`: With interaction

- `coe ~ poly(length, 2)`: Polynomial terms

### Getting results

Use [`collect()`](https://momx.github.io/Momstats/reference/collect.md)
to add predictions to your data:

    lm_fit <- boteft %>% stat_lm(coe ~ length)
    boteft_pred <- collect(lm_fit)  # Adds fitted values

Use
[`transduce()`](https://momx.github.io/Momstats/reference/transduce.md)
to predict at new positions:

    new_shapes <- transduce(lm_fit, tibble::tibble(length = seq(50, 150, by = 10)))

## See also

[`stats::lm()`](https://rdrr.io/r/stats/lm.html),
[`collect.stat_lm()`](https://momx.github.io/Momstats/reference/collect.stat_lm.md),
[`plot.stat_lm()`](https://momx.github.io/Momstats/reference/plot.stat_lm.md),
[`transduce.stat_lm()`](https://momx.github.io/Momstats/reference/transduce.stat_lm.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Shape regression
lm1 <- boteft %>% stat_lm(coe ~ length)
lm2 <- boteft %>% stat_lm(coe ~ length + type)
lm3 <- boteft %>% stat_lm(coe ~ length * type)

# Response regression
lm4 <- boteft %>% stat_lm(length ~ coe)

# Multivariate response
lm5 <- boteft %>% stat_lm(cbind(length, width) ~ coe)

# Plot results
plot(lm1)
plot(lm1, type = "residuals")

# Get predictions
boteft_pred <- collect(lm1)

# Transduce to new values
new_shapes <- transduce(lm1, tibble::tibble(length = seq(50, 150, by = 10)))
} # }
```
