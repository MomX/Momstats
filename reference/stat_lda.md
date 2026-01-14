# Linear Discriminant Analysis for morphometric data

Perform LDA on coefficient data from morphometric analyses. Supports
both bare response names (auto-detects all coe columns) and full formula
syntax.

## Usage

``` r
stat_lda(data, formula, ...)
```

## Arguments

- data:

  A tibble with coefficient columns

- formula:

  A formula or bare column name specifying response and predictors:

  - Bare name: `species` → auto-detects all coe columns

  - Formula: `species ~ .` → all coe columns

  - Formula: `species ~ . + length` → all coe + covariate

  - Formula: `species ~ coe1` → specific coe column

  - Formula: `species ~ coe1 + coe2 + length` → multiple specific

- ...:

  Additional arguments passed to
  [`MASS::lda()`](https://rdrr.io/pkg/MASS/man/lda.html)

## Value

An object of class `c("stat_lda", "momstats")` containing:

- `data`: Original tibble (unchanged)

- `model`: List containing both the LDA model and CV predictions:

  - All components from
    [`MASS::lda()`](https://rdrr.io/pkg/MASS/man/lda.html) (prior,
    counts, means, scaling, svd, etc.)

  - `cv_class`: Cross-validated predicted classes

  - `cv_posterior`: Cross-validated posterior probabilities

- `method`: "lda"

- `call`: The function call

- `formula`: Formula used

- `response_col`: Name of response column

- `coe_cols`: Names of coefficient columns used

- `covariate_cols`: Names of covariate columns used (if any)

- `predictor_cols`: All predictor column names

- `n_groups`: Number of groups

- `n_ld`: Number of discriminant functions

## Details

`stat_lda()` provides a unified interface for LDA on morphometric
coefficient data. It handles:

- **Flexible formula syntax**: Bare response names or full formulas

- **Automatic coefficient detection**: Finds all `coe` columns via class

- **Coefficient unfolding**: Expands list-columns to individual harmonic
  columns

- **Multiple coefficient types**: Can combine different coefficient
  columns

- **Covariates**: Add non-morphometric predictors (size, environmental
  variables)

- **Cross-validation**: Always runs with CV=TRUE for reliable
  predictions

### Formula syntax

The formula specifies response and predictors:

- `species` → Response only, auto-detect all coe columns

- `species ~ .` → All coe columns (explicit)

- `species ~ coe1` → Specific coe column

- `species ~ coe1 + coe2` → Multiple specific coe columns

- `species ~ . + length` → All coe plus a covariate

### Centering and scaling

LDA internally centers and scales the data as part of the algorithm. No
manual preprocessing is needed.

### Cross-validation

The function runs LDA twice:

1.  `CV = FALSE`: Builds the full model (discriminant functions,
    loadings)

2.  `CV = TRUE`: Generates cross-validated predictions

Both are stored in `$model` for easy access.

### Getting results

Use [`collect()`](https://momx.github.io/Momstats/reference/collect.md)
to add predictions to your data:

    lda <- boteft %>% stat_lda(type)
    boteft_pred <- collect(lda)  # Adds pred, prob columns

## See also

[`MASS::lda()`](https://rdrr.io/pkg/MASS/man/lda.html),
[`collect.stat_lda()`](https://momx.github.io/Momstats/reference/collect.stat_lda.md),
[`plot.stat_lda()`](https://momx.github.io/Momstats/reference/plot.stat_lda.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic LDA - auto-detects all coe columns
lda1 <- boteft %>% stat_lda(type)

# Explicit formula
lda2 <- boteft %>% stat_lda(type ~ .)

# Specific coe column
lda3 <- boteft %>% stat_lda(type ~ coe)

# With covariate
lda4 <- boteft %>% stat_lda(type ~ . + length)

# Add predictions to data
boteft_pred <- collect(lda1)  # Adds pred, prob

# Add LD scores too
boteft_ld <- collect(lda1, retain = 2)  # Adds pred, prob, LD1, LD2

# Plot results
plot(lda1)  # Discrimination plot (default)
plot(lda1, type = "loadings")
plot(lda1, color = type)  # Color by true groups
plot(lda1, labels = pred)  # Label by predictions
} # }
```
