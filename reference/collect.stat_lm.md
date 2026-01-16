# Collect fitted values from linear model

Extract fitted values and residuals from a linear model and add them to
a tibble.

## Usage

``` r
# S3 method for class 'stat_lm'
collect(x, data = NULL, fold = FALSE, residuals = FALSE, ...)
```

## Arguments

- x:

  A `stat_lm` object

- data:

  A tibble. If NULL, uses the original data from the model.

- fold:

  How to add fitted values:

  - `FALSE` (default): Add as columns matching response structure

  - `TRUE`: Fold into list-column (for coe responses)

  - Character: Name for folded column

- residuals:

  Logical. Should residuals be included? Default `FALSE`.

- ...:

  Additional arguments (reserved)

## Value

A tibble with fitted values (and optionally residuals) added

## Examples

``` r
if (FALSE) { # \dontrun{
lm_fit <- boteft %>% stat_lm(coe ~ length)

# Add fitted values
collect(lm_fit)

# With residuals
collect(lm_fit, residuals = TRUE)

# Fold into list-column
collect(lm_fit, fold = TRUE)
} # }
```
