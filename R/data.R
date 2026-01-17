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

#
# x <- Bacchus::Vitis_ref %>% Momocs::as_df()
# x %>% fold(VD=A1:D5, .class = c("eft", "coe")) %>% fold(VL=V21:V40, .class = c("eft", "coe")) %>% mutate(VL=map(VL, eft_name)) -> z
# class(z$VL) <- c("eft", "coe", "list")
# z <- z %>% mutate(ifp2=ifp %>% map_chr(~.x %>% str_split("_") %>% `[[`(1) %>% `[`(2))) %>%
#   mutate(accession=ifelse(status=="domesticated", as.character(variety), ifp2))
# z <- z %>% select(id, accession, status, length, usage=grow_usage, geo_group:VL)
#
# retain_these <- c(
# "Syrah",
# "Pinot blanc",
# "Pinot noir",
# "Merlot",
# "Grenache",
# "Cinsaut",
# "Cabernet-Sauvignon",
# "Clairette",
# "Mourvedre",
# "wPSL13",
# "wPSLH",
# "wLagMa1",
# "wCente3",
# "wGokce7")
#
# z %>% filter(accession %in% retain_these) %>%
#   group_by(accession) %>%
#   mutate(pip=1:n()) %>%
#   ungroup() %>%
#   select(id, pip, accession, pip, status, length, VD, VL) -> vitis
# usethis::use_data(vitis, overwrite = TRUE)


#' Grapevine pips outlines (two views) with elliptic Fourier transform coefficients
#'
#' A dataset containing 410 pip outlines (dorsal and lateral view) that have been preprocessed and
#' transformed using elliptic Fourier analysis (5 harmonics). This dataset is used for
#' demonstrating statistical methods in Momstats.
#'
#' @format A tibble with 410 rows and 5 variables:
#' * `id`: for the accession
#' * `pip`: Numeric identifying pip number
#' * `accession`: Character 14 names
#' * `status`: Factor with 2 levels (domesticated/wild)
#' * `length` : pip length measured before outline normalization
#' * `VD`: 5 eft harmonics for dorsal view
#' * `VL`: 5 eft harmonics for lateral view
#'
#' @source adapted from Bonhomme, Vincent; Ivorra, Sarah; Figueiral, Isabel; Pastor, Thierry; Terral, Jean-Frédéric; Bouby, Laurent (2021).
#' Vitis reference collection database (Sci Reports 2021). figshare. Dataset. https://doi.org/10.6084/m9.figshare.14170484.v1
#' @examples
#' # View structure
#' vitis
#'
"vitis"
