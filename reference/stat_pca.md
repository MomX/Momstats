# Principal Component Analysis for morphometric data

Perform PCA on coefficient data from morphometric analyses.

## Usage

``` r
stat_pca(data, formula = NULL, center = TRUE, scale = NULL, ...)
```

## Arguments

- data:

  A tibble with coefficient columns

- formula:

  A formula specifying predictors. Can be:

  - Missing: auto-detects single coe column

  - Bare column name: `coe`

  - Formula: `~ coe`, `~ coe + size`, `~ coe1 + coe2`

  - Use `coe` in formula to auto-detect coefficient columns

- center:

  Logical. Should data be centered? Default `TRUE`.

- scale:

  Logical or NULL. Should data be scaled to unit variance? If `NULL`
  (default), automatically determined based on predictor types.

- ...:

  Additional arguments passed to
  [`stats::prcomp()`](https://rdrr.io/r/stats/prcomp.html)

## Value

An object of class `c("stat_pca", "momstats")` containing:

- `data`: Original tibble (unchanged)

- `model`: The [`stats::prcomp()`](https://rdrr.io/r/stats/prcomp.html)
  object

- `method`: "pca"

- `call`: The function call

- `formula`: Formula used (if any)

- `coe_cols`: Names of coefficient columns used

- `covariate_cols`: Names of covariate columns used (if any)

- `predictor_cols`: All predictor column names

- `group_col`: NULL (no grouping for PCA)

- `variance_explained`: Proportion of variance explained by each PC

- `cumvar_explained`: Cumulative proportion of variance explained

- `center`: Logical, was centering applied

- `scale`: Logical, was scaling applied

## Details

`stat_pca()` provides a unified interface for PCA on morphometric
coefficient data. It handles:

- **Automatic coefficient detection**: Finds `coe` columns via class

- **Coefficient unfolding**: Expands list-columns to individual harmonic
  columns

- **Multiple coefficient types**: Can combine different coefficient
  columns

- **Covariates**: Add non-morphometric predictors (size, environmental
  variables)

- **Smart defaults**: Automatically determines centering and scaling

### Formula syntax

The formula specifies which predictors to use:

- `~ coe`: Use auto-detected coefficient column(s)

- `~ coe1 + coe2`: Use specific coefficient columns

- `~ coe + size`: Coefficient column plus a covariate

- Bare name or missing: auto-detect single coe column

### Centering and scaling

- `center = TRUE` (default): Variables are centered to mean = 0

- `scale = NULL` (default): Automatic determination:

  - Pure EFT coefficients: `scale = FALSE`

  - Mixed predictors or non-EFT: `scale = TRUE`

- `scale = TRUE/FALSE`: Override automatic behavior

### Getting results

Use [`collect()`](https://momx.github.io/Momstats/reference/collect.md)
to add PC scores to your data:

    pca <- boteft %>% stat_pca()
    boteft_with_pcs <- collect(pca, retain = 5)  # Add first 5 PCs

## See also

[`stats::prcomp()`](https://rdrr.io/r/stats/prcomp.html),
[`collect.stat_pca()`](https://momx.github.io/Momstats/reference/collect.stat_pca.md),
[`plot.stat_pca()`](https://momx.github.io/Momstats/reference/plot.stat_pca.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic PCA on coefficients
pca1 <- boteft %>% stat_pca()

# Explicit column
pca2 <- boteft %>% stat_pca(coe)

# Formula syntax
pca3 <- boteft %>% stat_pca(~ coe)

# With covariate
pca4 <- boteft %>% stat_pca(~ coe + length)

# Force scaling
pca5 <- boteft %>% stat_pca(~ coe, scale = TRUE)

# Add PC scores to data
boteft_pca <- collect(pca1, retain = 5)

# Or fold into single column
boteft_pca <- collect(pca1, retain = 5, fold = "pca_scores")

# Plot results
plot(pca1)  # Score plot (default)
plot(pca1, type = "scree")
plot(pca1, type = "loadings")
plot(pca1, color = type)  # Color by grouping variable
plot(pca1, labels = type)  # Text labels
plot(pca1, color = type, chull = TRUE)  # With convex hulls
} # }
```
