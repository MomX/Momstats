# Bottle outlines with elliptic Fourier transform coefficients

A dataset containing 40 bottle outlines that have been preprocessed and
transformed using elliptic Fourier analysis. This dataset is used for
demonstrating statistical methods in Momstats.

## Usage

``` r
boteft
```

## Format

A tibble with 40 rows and 5 variables:

- `coo`: List-column of class `out` containing outline coordinates (nx2
  matrices). Outlines have been centered, scaled to unit centroid size,
  aligned using PCA, and rotated so the rightmost point is at 0 radians.

- `type`: Factor with 2 levels ("whisky", "beer") indicating bottle
  type.

- `dummy`: Factor with 1 level ("a") - a dummy grouping variable for
  examples.

- `length`: Numeric vector containing the length of each outline along
  the major inertia axis (computed before standardization).

- `coe`: List-column of class `c("eft", "coe", "list")` containing
  elliptic Fourier transform coefficients. Each element is a named
  numeric vector with 24 coefficients (6 harmonics × 4 coefficients: A,
  B, C, D) representing the shape in Fourier space.

## Source

Obtained with:

    bot %>%
      dplyr::mutate(length = unlist(get_length(coo))) %>%
      coo_center() %>%
      coo_scale() %>%
      coo_align() %>%
      coo_slide_direction("right") %>%
      dplyr::mutate(coe = purrr::map(coo, eft, 6) %>%
                          as_coe() %>%
                          as_eft()) -> boteft

## Details

This dataset was created from the original `bot` dataset through the
following preprocessing pipeline:

1.  **Length measurement**: Original outline lengths were computed and
    stored

2.  **Centering**: Outlines translated to origin (centroid at 0,0)

3.  **Scaling**: Outlines scaled to unit centroid size

4.  **Alignment**: Outlines aligned along major inertia axis using PCA

5.  **Rotation**: Outlines rotated so rightmost point is at 0 radians

6.  **Fourier transform**: Elliptic Fourier analysis with 6 harmonics

The preprocessing standardizes size, position, and orientation, making
the outlines directly comparable for statistical analysis. The elliptic
Fourier coefficients provide a compact representation of shape that can
be used with multivariate statistical methods like PCA, LDA, and
clustering.

### Fourier Coefficients

Each coefficient vector contains 24 values organized as:

- A1-A6: First harmonic coefficients (x-axis cosine terms)

- B1-B6: Second harmonic coefficients (x-axis sine terms)

- C1-C6: Third harmonic coefficients (y-axis cosine terms)

- D1-D6: Fourth harmonic coefficients (y-axis sine terms)

Higher harmonics capture finer shape details, while lower harmonics
represent the overall shape structure.

## See also

- [Momocs2::bot](https://momx.github.io/Momocs2/reference/bot.html) for
  the original unprocessed bottle outlines

- [`Momocs2::unfold()`](https://momx.github.io/Momocs2/reference/unfold.html)
  to expand coefficient columns for analysis

- [`Momocs2::fold()`](https://momx.github.io/Momocs2/reference/fold.html)
  for the opposite transformation

## Examples

``` r
# View structure
boteft
#> # A tibble: 40 × 5
#>    coo       type   dummy length coe  
#>    <out>     <fct>  <fct>  <dbl> <eft>
#>  1 (138 x 2) whisky a      1088. <NA> 
#>  2 (168 x 2) whisky a       994. <NA> 
#>  3 (189 x 2) whisky a       644. <NA> 
#>  4 (129 x 2) whisky a       806. <NA> 
#>  5 (152 x 2) whisky a       886. <NA> 
#>  6 (161 x 2) whisky a       606. <NA> 
#>  7 (124 x 2) whisky a       865. <NA> 
#>  8 (126 x 2) whisky a       765. <NA> 
#>  9 (183 x 2) whisky a       742. <NA> 
#> 10 (193 x 2) whisky a      1048. <NA> 
#> # ℹ 30 more rows

# Examine coefficient structure
class(boteft$coe)
#> [1] "eft"  "coe"  "list"
class(boteft$coe[[1]])
#> [1] "eft"     "numeric"
boteft$coe[[1]]
#>            A1            A2            A3            A4            A5 
#>  1.3427539129  0.0090863444  0.1248632443  0.0183772998  0.0314900219 
#>            A6            B1            B2            B3            B4 
#>  0.0113182436  0.0469665913  0.0003981801  0.0136025769  0.0019011705 
#>            B5            B6            C1            C2            C3 
#>  0.0057056253  0.0022844660  0.0134712409 -0.0021381830  0.0126151941 
#>            C4            C5            C6            D1            D2 
#> -0.0139137033  0.0110247077  0.0084896655 -0.3943943267  0.0618459428 
#>            D3            D4            D5            D6 
#> -0.0694671431  0.0465708687 -0.0526203987 -0.0140823150 
#> attr(,"class")
#> [1] "eft"     "numeric"

# Unfold coefficients for statistical analysis
boteft_unfolded <- unfold(boteft, coe)

# Or without prefix
boteft_unfolded <- unfold(boteft, coe, .prefix = "")

# Principal component analysis on coefficients
# pca_result <- stat_pca(boteft)

# Linear discriminant analysis by bottle type
# lda_result <- stat_lda(boteft, type)
```
