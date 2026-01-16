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

- `coo`: List-column of class `c("out", "coo", "list")` containing
  outline coordinates (nx2 matrices). Outlines have been centered,
  scaled to unit centroid size, aligned using PCA, and rotated so the
  rightmost point is at 0 radians.

- `coe`: List-column of class `c("eft", "coe", "list")` containing
  elliptic Fourier transform coefficients. Each element is a named
  numeric vector with 24 coefficients (6 harmonics × 4 coefficients: A,
  B, C, D) representing the shape in Fourier space.

- `type`: Factor with 2 levels ("whisky", "beer") indicating bottle
  type.

- `fake`: Factor with 1 level ("a") - a dummy grouping variable for
  examples.

- `price` : Numeric with random prices

- `length`: Numeric vector containing the length of each outline along
  the major inertia axis (computed before standardization, also random
  because it ultimately depends on the size of original images in
  pixels...)

## Source

Obtained with:

    bot %>%
      mutate(length=unlist(get_length(coo))) %>%
      coo_center() %>% coo_scale() %>% coo_align() %>%
      coo_slide_direction("right") %>%
      eft(6) %>% relocate(coe, .after=coo) -> boteft

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
#> # A tibble: 40 × 7
#>    id           coo       coe   type   fake  price length
#>    <chr>        <out>     <eft> <fct>  <fct> <dbl>  <dbl>
#>  1 brahma       (138 x 2) <24>  whisky a       3    1088.
#>  2 caney        (168 x 2) <24>  whisky a       1.2   994.
#>  3 chimay       (189 x 2) <24>  whisky a       3.8   644.
#>  4 corona       (129 x 2) <24>  whisky a       2.6   806.
#>  5 deusventrue  (152 x 2) <24>  whisky a       1.1   886.
#>  6 duvel        (161 x 2) <24>  whisky a       3.1   606.
#>  7 franziskaner (124 x 2) <24>  whisky a       2.6   865.
#>  8 grimbergen   (126 x 2) <24>  whisky a       2.9   765.
#>  9 guiness      (183 x 2) <24>  whisky a       1.2   742.
#> 10 hoegardeen   (193 x 2) <24>  whisky a       3.6  1048.
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
