# Predict method for PCA

Project new data onto principal components from a fitted PCA model.

## Usage

``` r
# S3 method for class 'stat_pca'
predict(object, newdata, retain = NULL, fold = FALSE, .collect = TRUE, ...)
```

## Arguments

- object:

  A `stat_pca` object

- newdata:

  A tibble with the same predictor columns as training data

- retain:

  How many PCs to return:

  - `NULL` (default): All PCs

  - Integer (e.g., `5`): First N PCs

  - Numeric 0-1 (e.g., `0.95`): PCs explaining this proportion of
    variance

- fold:

  How to return PC scores:

  - `FALSE` (default): Add as separate columns (`PC1`, `PC2`, ...)

  - `TRUE`: Fold into single list-column named `"pca"`

  - Character: Fold into single list-column with this name

- .collect:

  Logical. Should predictions be added to `newdata` (TRUE, default) or
  returned as a standalone tibble (FALSE)?

- ...:

  Additional arguments (reserved)

## Value

If `.collect = TRUE`, returns `newdata` with PC scores added. If
`.collect = FALSE`, returns a tibble with identifier columns (non-coe)
and PC scores only.

## Examples

``` r
if (FALSE) { # \dontrun{
# Train PCA
pca <- boteft %>% stat_pca()

# Predict on new data
new_scores <- predict(pca, new_data)

# Pipe the model
new_scores <- training %>%
  stat_pca() %>%
  predict(testing)
} # }
```
