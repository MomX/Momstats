# Collect predictions and LD scores from LDA results

Extract predictions, posterior probabilities, and optionally LD scores
from an LDA result and add them to a tibble.

## Usage

``` r
# S3 method for class 'stat_lda'
collect(x, data = NULL, retain = FALSE, fold = FALSE, ...)
```

## Arguments

- x:

  A `stat_lda` object

- data:

  A tibble. If NULL, uses the original data from the LDA.

- retain:

  How many LDs to retain:

  - `FALSE` (default): No LD scores, only predictions

  - `TRUE`: All LDs

  - Integer (e.g., `2`): First N LDs

  - Numeric 0-1 (e.g., `0.95`): LDs explaining this proportion of
    variance

- fold:

  How to add LD scores (only if retain != FALSE):

  - `FALSE` (default): Add as separate columns (`LD1`, `LD2`, ...)

  - `TRUE`: Fold into single list-column named `"lda"`

  - Character: Fold into single list-column with this name

- ...:

  Additional arguments (reserved)

## Value

A tibble with predictions added:

- `pred`: Predicted class (factor)

- `prob`: Posterior probability of predicted class (numeric)

- `LD1`, `LD2`, ... : LD scores (if `retain` is not FALSE)

## Details

The function always adds cross-validated predictions (`pred`) and the
posterior probability of the predicted class (`prob`). These come from
leave-one-out cross-validation performed during model fitting.

Optionally, discriminant scores can be added via the `retain` parameter:

- `retain = FALSE`: No LD scores (default, fastest)

- `retain = TRUE`: All LD scores

- `retain = 2`: First 2 LDs

- `retain = 0.95`: LDs explaining 95% of between-group variance

When `fold = TRUE` or a character name, LD scores are stored as a
list-column with class `c("lda", "coe")`, making them usable in
downstream analyses.

## Examples

``` r
if (FALSE) { # \dontrun{
lda <- boteft %>% stat_lda(type)

# Add predictions only (default)
collect(lda)

# Add predictions + all LD scores
collect(lda, retain = TRUE)

# Add predictions + first 2 LDs
collect(lda, retain = 2)

# Add predictions + LDs explaining 95% variance
collect(lda, retain = 0.95)

# Fold LD scores into list-column
collect(lda, retain = 2, fold = TRUE)
collect(lda, retain = 2, fold = "lda_scores")
} # }
```
