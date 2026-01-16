# Predict method for LDA

Predict class membership and discriminant scores for new data using a
fitted LDA model.

## Usage

``` r
# S3 method for class 'stat_lda'
predict(object, newdata, retain = FALSE, fold = FALSE, .collect = TRUE, ...)
```

## Arguments

- object:

  A `stat_lda` object

- newdata:

  A tibble with the same predictor columns as training data

- retain:

  How many LDs to return:

  - `FALSE` (default): No LD scores, only predictions

  - `TRUE`: All LDs

  - Integer (e.g., `2`): First N LDs

  - Numeric 0-1 (e.g., `0.95`): LDs explaining this proportion of
    variance

- fold:

  How to return LD scores (only if retain != FALSE):

  - `FALSE` (default): Add as separate columns (`LD1`, `LD2`, ...)

  - `TRUE`: Fold into single list-column named `"lda"`

  - Character: Fold into single list-column with this name

- .collect:

  Logical. Should predictions be added to `newdata` (TRUE, default) or
  returned as a standalone tibble (FALSE)?

- ...:

  Additional arguments (reserved)

## Value

If `.collect = TRUE`, returns `newdata` with predictions added. If
`.collect = FALSE`, returns a tibble with identifier columns
(non-predictors) and predictions only. Always includes:

- `pred`: Predicted class (factor)

- `prob`: Posterior probability of predicted class (numeric)

- `LD1`, `LD2`, ... : LD scores (if `retain` is not FALSE)

## Examples

``` r
if (FALSE) { # \dontrun{
# Train LDA
lda <- boteft %>% stat_lda(type)

# Predict on new data
new_preds <- predict(lda, new_data)

# Pipe the model
new_preds <- training %>%
  stat_lda(type) %>%
  predict(testing)
} # }
```
