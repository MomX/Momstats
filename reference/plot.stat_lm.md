# Plot linear model results

Visualize linear model diagnostics or predictions.

## Usage

``` r
# S3 method for class 'stat_lm'
plot(x, type = c("residuals", "qq", "fitted", "predictions"), which = 1, ...)
```

## Arguments

- x:

  A `stat_lm` object

- type:

  Character. Type of plot:

  - `"residuals"`: Residual plot

  - `"qq"`: Q-Q plot

  - `"fitted"`: Fitted vs observed

  - `"predictions"`: Predictions along predictor (shape regression only)

- which:

  For diagnostic plots, which response/harmonic to plot (integer)

- ...:

  Additional arguments passed to plot functions

## Value

NULL (invisibly). Draws plot as side effect.

## Examples

``` r
if (FALSE) { # \dontrun{
lm_fit <- boteft %>% stat_lm(coe ~ length)

plot(lm_fit, type = "residuals")
plot(lm_fit, type = "qq")
plot(lm_fit, type = "predictions")
} # }
```
