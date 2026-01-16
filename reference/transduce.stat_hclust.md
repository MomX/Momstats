# Transduce hierarchical clustering nodes to shapes

Reconstruct shapes at cluster centers or at any node in the tree.

## Usage

``` r
# S3 method for class 'stat_hclust'
transduce(object, positions)
```

## Arguments

- object:

  A `stat_hclust` object

- positions:

  A tibble specifying positions. Can contain:

  - `cluster` column: cluster numbers (if tree was cut)

  - `node` column: node numbers in the tree (tips are 1:n, internal
    nodes are n+1:...)

## Value

A tibble with reconstructed coefficients and shapes

## Details

For hierarchical clustering, transduce can reconstruct shapes at:

**Cluster centers** (if tree was cut with k or h):

    transduce(hc, tibble(cluster = 1:4))

Returns the mean shape of each cluster.

**Any node in the tree**:

    transduce(hc, tibble(node = c(45, 50, 55)))

Returns the mean shape of all descendants of that node. This allows
visualization of "ancestral" shapes at internal nodes.

Node numbering follows
[`stats::hclust()`](https://rdrr.io/r/stats/hclust.html) convention:

- Nodes 1 to n: original observations (tips)

- Nodes n+1 to 2n-1: internal nodes from merges

## Examples

``` r
if (FALSE) { # \dontrun{
hc <- boteft %>% stat_hclust(k = 4)

# Cluster centers
centers <- transduce(hc, tibble(cluster = 1:4))

# Shapes at internal nodes
hc2 <- boteft %>% stat_hclust()
internal <- transduce(hc2, tibble(node = c(45, 50, 55, 60)))

# Root node (all data)
root <- transduce(hc2, tibble(node = 2 * nrow(boteft) - 1))
} # }
```
