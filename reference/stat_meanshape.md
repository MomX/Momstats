# Calculate mean (or aggregated) shapes

Compute mean or other aggregated shapes, either globally or by groups.

## Usage

``` r
stat_meanshape(data, formula = NULL, .fn = mean, ...)
```

## Arguments

- data:

  A tibble with coefficient columns

- formula:

  A formula or bare column name specifying grouping:

  - Missing/NULL: Global mean across all observations

  - Bare name: `status` → group by this column, use all coe columns

  - Formula: `status ~ .` → explicit all coe columns

  - Formula: `status ~ VD` → specific coe column(s)

  - Formula: `~ status` → one-sided formula (group by status)

- .fn:

  Aggregation function. Default `mean`. Can be `median` or any function
  that takes a matrix column and returns a scalar.

- ...:

  Additional arguments (reserved for future use)

## Value

An object of class `c("stat_meanshape", "momstats")` containing:

- `data`: Original tibble (unchanged)

- `summary`: Tibble with group(s), n, and mean coe column(s)

- `method`: "meanshape"

- `call`: The function call

- `formula`: Formula used (if any)

- `group_cols`: Grouping column name(s) (NULL if global)

- `coe_cols`: Coefficient column names

- `.fn`: Aggregation function name

## Details

`stat_meanshape()` provides a unified interface for aggregating shapes.
It handles:

- **Global aggregation**: No grouping, single mean shape

- **Grouped aggregation**: By one or more grouping variables

- **Flexible syntax**: Bare names or formula syntax

- **Multiple coe columns**: Aggregates all coe columns simultaneously

- **Custom aggregation**: mean (default), median, or custom function

- **Automatic inverse**: Computes reconstructed shapes via
  `coe_class_i()`

### Formula syntax

The formula specifies grouping and optionally which coe columns:

- Missing: Global mean of all coe columns

- `status`: Group by status, all coe columns

- `~ status`: Same as above

- `status ~ .`: Explicit all coe columns

- `status ~ VD`: Only VD coe column

- `status ~ VD + VL`: Multiple specific coe columns

### Aggregation functions

- `.fn = mean` (default): Arithmetic mean

- `.fn = median`: Median shape

- Custom: Any function taking a numeric vector and returning a scalar

### Getting results

Use [`collect()`](https://momx.github.io/Momstats/reference/collect.md)
to add mean shapes back to original data (matched by group):

    ms <- vitis %>% stat_meanshape(status)
    vitis_with_means <- collect(ms)  # Adds mean_VD, mean_VL, etc.

## See also

[`collect.stat_meanshape()`](https://momx.github.io/Momstats/reference/collect.stat_meanshape.md),
[`plot.stat_meanshape()`](https://momx.github.io/Momstats/reference/plot.stat_meanshape.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Global mean
ms1 <- vitis %>% stat_meanshape()

# By group
ms2 <- vitis %>% stat_meanshape(status)

# Formula syntax
ms3 <- vitis %>% stat_meanshape(~ status)
ms4 <- vitis %>% stat_meanshape(status ~ .)

# Specific coe column
ms5 <- vitis %>% stat_meanshape(status ~ VD)

# Median instead of mean
ms6 <- vitis %>% stat_meanshape(status, .fn = median)

# Add mean shapes to data
vitis_means <- collect(ms2)

# Plot mean shapes
plot(ms2)  # Mean shapes only
plot(ms2, with_data = TRUE)  # Mean shapes + contributing data
} # }
```
