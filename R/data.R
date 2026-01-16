# # on fresh session:
# library(Momocs2)
# library(tidyverse)
#
# bot %>%
#   mutate(length=unlist(get_length(coo))) %>%
#   coo_center() %>% coo_scale() %>% coo_align() %>%
#   coo_slide_direction("right") %>%
#   eft(6) %>% relocate(coe, .after=coo) -> boteft
# usethis::use_data(boteft, overwrite = TRUE)

#' Bottle outlines with elliptic Fourier transform coefficients
#'
#' A dataset containing 40 bottle outlines that have been preprocessed and
#' transformed using elliptic Fourier analysis. This dataset is used for
#' demonstrating statistical methods in Momstats.
#'
#' @format A tibble with 40 rows and 5 variables:
#' * `coo`: List-column of class `c("out", "coo", "list")` containing outline coordinates
#'   (nx2 matrices). Outlines have been centered, scaled to unit centroid size,
#'   aligned using PCA, and rotated so the rightmost point is at 0 radians.
#' * `coe`: List-column of class `c("eft", "coe", "list")` containing
#'   elliptic Fourier transform coefficients. Each element is a named numeric
#'   vector with 24 coefficients (6 harmonics × 4 coefficients: A, B, C, D)
#'   representing the shape in Fourier space.
#' * `type`: Factor with 2 levels ("whisky", "beer") indicating bottle type.
#' * `fake`: Factor with 1 level ("a") - a dummy grouping variable for examples.
#' * `price` : Numeric with random prices
#' * `length`: Numeric vector containing the length of each outline along the
#'   major inertia axis (computed before standardization, also random because it
#'   ultimately depends on the size of original images in pixels...)
#'
#' @source Obtained with:
#' ```r
#' bot %>%
#'   mutate(length=unlist(get_length(coo))) %>%
#'   coo_center() %>% coo_scale() %>% coo_align() %>%
#'   coo_slide_direction("right") %>%
#'   eft(6) %>% relocate(coe, .after=coo) -> boteft
#' ```
#'
#' @seealso
#' * [Momocs2::bot] for the original unprocessed bottle outlines
#' * [Momocs2::unfold()] to expand coefficient columns for analysis
#' * [Momocs2::fold()] for the opposite transformation
# * [stat_pca()], [stat_lda()] for statistical analysis of coefficients
#'
#' @examples
#' # View structure
#' boteft
#'
#' # Examine coefficient structure
#' class(boteft$coe)
#' class(boteft$coe[[1]])
#' boteft$coe[[1]]
#'
#' # Unfold coefficients for statistical analysis
#' boteft_unfolded <- unfold(boteft, coe)
#'
#' # Or without prefix
#' boteft_unfolded <- unfold(boteft, coe, .prefix = "")
#'
#' # Principal component analysis on coefficients
#' # pca_result <- stat_pca(boteft)
#'
#' # Linear discriminant analysis by bottle type
#' # lda_result <- stat_lda(boteft, type)
#'
"boteft"
