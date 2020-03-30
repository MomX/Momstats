# outlines ------------------------------------------------
#' Elliptical Fourier transforms coefficients from bottles dataset
#'
#' Obtained with `bot2 %>% efourier(6) %>% dplyr::mutate(q=rnorm(40, 5)) %>% dplyr::select(coe, type, fake, q)`
#' @docType data
#' @name bot_eft
#' @family outlines
#' @format A `coe_tbl` with:
#'  * `coe` (`coe_list`): eft coefficients
#'  * `type` (`factor`): bottle type
#'  * `fake` (`factor`): bottle type
#'
#' @source Borrowed default shapes from (c) Adobe Photoshop. Please do not send me to jail.
#' @examples
#' bot_eft
NULL

# bot2 %>% efourier(6) %>%
#   dplyr::mutate(q=rnorm(40, 5)) dplyr::select(coe, type, fake, q) %>%
#   usethis::use_data(bot_eft, overwrite = TRUE)

