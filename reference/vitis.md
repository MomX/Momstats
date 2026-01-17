# Grapevine pips outlines (two views) with elliptic Fourier transform coefficients

A dataset containing 410 pip outlines (dorsal and lateral view) that
have been preprocessed and transformed using elliptic Fourier analysis
(5 harmonics). This dataset is used for demonstrating statistical
methods in Momstats.

## Usage

``` r
vitis
```

## Format

A tibble with 410 rows and 5 variables:

- `id`: for the accession

- `pip`: Numeric identifying pip number

- `accession`: Character 14 names

- `status`: Factor with 2 levels (domesticated/wild)

- `length` : pip length measured before outline normalization

- `VD`: 5 eft harmonics for dorsal view

- `VL`: 5 eft harmonics for lateral view

## Source

adapted from Bonhomme, Vincent; Ivorra, Sarah; Figueiral, Isabel;
Pastor, Thierry; Terral, Jean-Frédéric; Bouby, Laurent (2021). Vitis
reference collection database (Sci Reports 2021). figshare. Dataset.
https://doi.org/10.6084/m9.figshare.14170484.v1

## Examples

``` r
# View structure
vitis
#> # A tibble: 420 × 7
#>    id      pip accession          status       length VD    VL   
#>    <chr> <int> <chr>              <fct>         <dbl> <eft> <eft>
#>  1 0036      1 Cabernet-Sauvignon domesticated   5.78 <20>  <20> 
#>  2 0036      2 Cabernet-Sauvignon domesticated   5.56 <20>  <20> 
#>  3 0036      3 Cabernet-Sauvignon domesticated   5.40 <20>  <20> 
#>  4 0036      4 Cabernet-Sauvignon domesticated   5.75 <20>  <20> 
#>  5 0036      5 Cabernet-Sauvignon domesticated   5.83 <20>  <20> 
#>  6 0036      6 Cabernet-Sauvignon domesticated   5.75 <20>  <20> 
#>  7 0036      7 Cabernet-Sauvignon domesticated   5.54 <20>  <20> 
#>  8 0036      8 Cabernet-Sauvignon domesticated   5.66 <20>  <20> 
#>  9 0036      9 Cabernet-Sauvignon domesticated   5.55 <20>  <20> 
#> 10 0036     10 Cabernet-Sauvignon domesticated   5.42 <20>  <20> 
#> # ℹ 410 more rows
```
