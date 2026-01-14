# Collect PC scores from PCA results

Extract PC scores from a PCA result and add them to a tibble.

## Usage

``` r
# S3 method for class 'stat_pca'
collect(x, data = NULL, retain = NULL, fold = FALSE, ...)
```

## Arguments

- x:

  A `stat_pca` object

- data:

  A tibble. If NULL, uses the original data from the PCA.

- retain:

  How many PCs to retain:

  - `NULL` (default): All PCs

  - Integer (e.g., `5`): First N PCs

  - Numeric 0-1 (e.g., `0.95`): PCs explaining this proportion of
    variance

- fold:

  How to add PC scores:

  - `FALSE` (default): Add as separate columns (`PC1`, `PC2`, ...)

  - `TRUE`: Fold into single list-column named `"pca"`

  - Character: Fold into single list-column with this name

- ...:

  Additional arguments (reserved)

## Value

A tibble with PC scores added

## Details

The `retain` parameter allows flexible selection of PCs:

- `retain = NULL`: Keep all PCs

- `retain = 5`: Keep first 5 PCs

- `retain = 0.95`: Keep PCs explaining 95% of variance

- `retain = 1`: Keep only PC1

When `fold = TRUE` or a character name, PC scores are stored as a
list-column with class `c("pca", "coe")`, making them usable in
downstream analyses.

## Examples

``` r
if (FALSE) { # \dontrun{
pca <- boteft %>% stat_pca()

# Add all PCs
collect(pca)

# Add first 5 PCs
collect(pca, retain = 5)

# Add PCs explaining 95% variance
collect(pca, retain = 0.95)

# Fold into list-column
collect(pca, retain = 5, fold = TRUE)
collect(pca, retain = 5, fold = "pca_scores")
} # }
```
