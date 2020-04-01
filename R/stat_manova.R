# stat_manova =============================================

# utils ---------------------------------------------------
.cross <- function(x){
  if (is.factor(x))
    x <- levels(x)
  x %>% combn(2) %>% t %>% tibble::as_tibble() %>% `colnames<-`(c("pair1", "pair2"))
}

#
# set.seed(2329)
# k=100
#
# df <- tibble::tibble(
#   # q variables, some collinear,
#   a=rnorm(k), b=jitter(a), c=runif(k),
#   # now constants, collinear or not
#   d=pi, e=jitter(d), f=3,
#   # factors now
#   foo1=factor("a"),
#   foo2=sample(LETTERS[1:12], size=k, replace = TRUE) %>% factor(),
#   foo3=sample(c("yes", "no"), size=k, replace=TRUE) %>% factor(),
# )
# e_NA <- df$e
# e_NA[c(5, 12)] <- NA
#
# foo2_NA <- df$foo2
# foo2_NA[c(7, 34)] <- NA
#
# df <- df %>% dplyr::mutate(e_NA=e_NA, foo2_NA=foo2_NA)
#
# x <- stat_lda_prepare(df, foo2_NA, a:f)
# df <- x$df

# stat_manova ---------------------------------------------
#' Multivariate analysis of variance
#'
#' Calculations are delegated to [stats::manova]
#'
#' @param x [tibble][tibble::tibble-package]
#' @param f `factor` used as right hand side in [stats::manova]'s formula
#' @param ... columns used as left hand side in [stats::manova]'s formula
#' @param test forwarded to [broom::tidy.manova]'s `test`. One of "Pillai" (Pillai's trace, default),
#' "Wilks" (Wilk's lambda), "Hotelling-Lawley" (Hotelling-Lawley trace) or
#' "Roy" (Roy's greatest root) indicating which test statistic should be used.
#'
#' @return see [broom::tidy.manova]'s Value section
#' @family manova
#' @export
stat_manova <- function(x, f, ..., test = "Pillai"){
  # tidyeval
  f_enquo   <- enquo(f)  ### todo: handles for formula here ??
  coe_enquo <- rlang::expr(c(...))
  coe_names <- tidyselect::eval_select(coe_enquo, x) ### ... should add to here

  # f first, coe then, drop others
  # f will be handled positionnally after (1st col)
  # so that it's simple and does not require colname handling
  coe <- dplyr::select(x, !!coe_enquo) %>% as.matrix()
  f <- dplyr::pull(x, !!f_enquo)
  stats::manova(coe ~ f) %>% broom::tidy(test=test)
}

# stat_manova(df, foo2_NA, a:e)

#' @export
#' @describeIn stat_manova pairwise manova for every pair of levels
stat_manova_pairwise <- function(x, f, ..., test = "Pillai"){
  # tidyeval
  f_enquo   <- enquo(f)  ### todo: handles for formula here ??
  coe_enquo <- rlang::expr(c(...))
  coe_names <- tidyselect::eval_select(coe_enquo, x) ### ... should add to here


  f <- dplyr::pull(x, !!f_enquo) %>% forcats::fct_drop()
  f_df <- .cross(f)
  n <- tidy_res <- vector("list", nrow(f_df))

  # loop over all combinations
  for (i in 1:nrow(f_df)){
    levels_i      <- unlist(f_df[i, ])
    x_i           <- dplyr::filter(x, !!f_enquo %in% levels_i)
    coe_i         <- dplyr::select(x_i, !!coe_enquo) %>% as.matrix()
    f_i           <- dplyr::pull(x_i, !!f_enquo) %>% forcats::fct_drop()
    n[[i]]        <- purrr::map_dbl(levels_i, ~sum(f_i==.x)) %>% `names<-`(c("pair1_n", "pair2_n"))
    tidy_res[[i]] <- stats::manova(coe_i ~ f_i) %>%
      broom::tidy(test=test) %>% dplyr::slice(1) %>% dplyr::select(-term)
  }

  # prepare and return this beauty
  dplyr::bind_cols(f_df,
                   tibble::tibble(pair1_n=map_dbl(n, 1), pair2_n=map_dbl(n, 2)),
                   dplyr::bind_rows(tidy_res))
}

# df %>%
#   dplyr::filter(foo2_NA %in% LETTERS[1:3]) %>%
#   stat_manova_pairwise(foo2_NA, a:e, test="Wilks")
#


