# bot %>% dplyr::mutate(length=unlist(get_length(coo))) %>% coo_center %>% coo_scale %>% coo_align() %>% coo_slide_direction("right") %>%
#   dplyr::mutate(coe=purrr:::map(coo, eft, 6) %>% as_coe() %>% as_eft()) -> boteft
# usethis::use_data(boteft)

#' Bottle outlines with elliptic Fourier transform coefficients
#'
#' A dataset containing 40 bottle outlines that have been preprocessed and
#' transformed using elliptic Fourier analysis. This dataset is used for
#' demonstrating statistical methods in Momstats.
#'
#' @format A tibble with 40 rows and 5 variables:
#' * `coo`: List-column of class `out` containing outline coordinates
#'   (nx2 matrices). Outlines have been centered, scaled to unit centroid size,
#'   aligned using PCA, and rotated so the rightmost point is at 0 radians.
#' * `type`: Factor with 2 levels ("whisky", "beer") indicating bottle type.
#' * `dummy`: Factor with 1 level ("a") - a dummy grouping variable for examples.
#' * `length`: Numeric vector containing the length of each outline along the
#'   major inertia axis (computed before standardization).
#' * `coe`: List-column of class `c("eft", "coe", "list")` containing
#'   elliptic Fourier transform coefficients. Each element is a named numeric
#'   vector with 24 coefficients (6 harmonics × 4 coefficients: A, B, C, D)
#'   representing the shape in Fourier space.
#'
#' @details
#' This dataset was created from the original `bot` dataset through the
#' following preprocessing pipeline:
#'
#' 1. **Length measurement**: Original outline lengths were computed and stored
#' 2. **Centering**: Outlines translated to origin (centroid at 0,0)
#' 3. **Scaling**: Outlines scaled to unit centroid size
#' 4. **Alignment**: Outlines aligned along major inertia axis using PCA
#' 5. **Rotation**: Outlines rotated so rightmost point is at 0 radians
#' 6. **Fourier transform**: Elliptic Fourier analysis with 6 harmonics
#'
#' The preprocessing standardizes size, position, and orientation, making the
#' outlines directly comparable for statistical analysis. The elliptic Fourier
#' coefficients provide a compact representation of shape that can be used with
#' multivariate statistical methods like PCA, LDA, and clustering.
#'
#' ## Fourier Coefficients
#'
#' Each coefficient vector contains 24 values organized as:
#' * A1-A6: First harmonic coefficients (x-axis cosine terms)
#' * B1-B6: Second harmonic coefficients (x-axis sine terms)
#' * C1-C6: Third harmonic coefficients (y-axis cosine terms)
#' * D1-D6: Fourth harmonic coefficients (y-axis sine terms)
#'
#' Higher harmonics capture finer shape details, while lower harmonics represent
#' the overall shape structure.
#'
#' @source Obtained with:
#' ```r
#' bot %>%
#'   dplyr::mutate(length = unlist(get_length(coo))) %>%
#'   coo_center() %>%
#'   coo_scale() %>%
#'   coo_align() %>%
#'   coo_slide_direction("right") %>%
#'   dplyr::mutate(coe = purrr::map(coo, eft, 6) %>%
#'                       as_coe() %>%
#'                       as_eft()) -> boteft
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
"boteft"
