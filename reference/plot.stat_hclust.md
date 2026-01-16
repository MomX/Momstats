# Plot hierarchical clustering tree

Visualize hierarchical clustering with various tree layouts.

## Usage

``` r
# S3 method for class 'stat_hclust'
plot(
  x,
  type = c("unrooted", "dendrogram", "fan", "phylogram"),
  color = NULL,
  tip_labels = NULL,
  ...
)
```

## Arguments

- x:

  A `stat_hclust` object

- type:

  Character. Type of plot:

  - `"unrooted"`: Unrooted phylogram (default, requires ape)

  - `"dendrogram"`: Classic dendrogram

  - `"fan"`: Fan/circular tree (requires ape)

  - `"phylogram"`: Rectangular phylogram (requires ape)

- color:

  Column name (bare or quoted) for coloring tips. If NULL and tree was
  cut, colors by cluster assignment.

- tip_labels:

  Logical. Show tip labels? Default TRUE for small trees (\<50 tips).

- ...:

  Additional arguments passed to plotting functions

## Value

NULL (invisibly). Draws plot as side effect.

## Details

Most plot types require the `ape` package. If not installed, falls back
to base dendrogram plot.

## Examples

``` r
if (FALSE) { # \dontrun{
hc <- boteft %>% stat_hclust(k = 4)

# Unrooted tree (default)
plot(hc)

# Color by original grouping
plot(hc, color = type)

# Color by clusters
plot(hc)  # Automatic if tree was cut

# Classic dendrogram
plot(hc, type = "dendrogram")

# Fan tree
plot(hc, type = "fan")
} # }
```
