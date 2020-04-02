# utils ---------------------------------------------------

# internal to turn matrics to CV_tbl
.CV_tbl <- function(x){
  # in case unnammed matrices are passed
  x <- as.table(x)
  names(dimnames(x)) <- c("actual", "predicted")
  x %>%
    tibble::as_tibble() %>%
    dplyr::mutate(n_total=sum(.data$n)) %>%
    # class acc
    dplyr::group_by(.data$actual) %>%
    dplyr::mutate(n_class=sum(.data$n),
                  prop_class=.data$n_class/.data$n_total) %>%
    dplyr::group_by(.data$actual, .data$predicted) %>%
    # case acc
    dplyr::mutate(prop=.data$n/.data$n_class) %>%
    dplyr::ungroup() %>%
    # cosmetics
    dplyr::select(.data$actual, .data$predicted,
                  .data$n, .data$prop,
                  .data$n_class, .data$prop_class, .data$n_total)
}

# stat_lda_prepare ----------------------------------------

#' Linear discriminant analysis preparation
#'
#' Removes NAs, constant and collinear columns.
#' Useful to speed up calculations when used with [stat_lda0]
#'
#' @param x [tibble][tibble::tibble-package], typically a [coe_tbl]
#' @param f column specifying grouping factor
#' @param ... columns specifying those to use
#' @param constant_below_var `numeric` threshold under which columns
#' will be considered constant (using [stats::var]) and dropped (default to `1e-5`)
#' @param collinear_above_cor `numeric` threshold above which columns
#' will be considered collinear (using [stats::cor]), and dropped (default to `1-1e-5`)
#'
#' @family lda
#' @export
stat_lda_prepare <- function(x,
                             f,
                             ...,
                             constant_below_var  = 1e-5,
                             collinear_above_cor = 1-1e-5){
  # tidyeval
  f_enquo   <- enquo(f)  ### todo: handles for formula here ??
  coe_enquo <- rlang::expr(c(...))
  coe_names <- tidyselect::eval_select(coe_enquo, x) ### ... should add to here

  # f first, coe then, drop others
  # f will be handled positionnally after (1st col)
  # so that it's simple and does not require colname handling
  df <- dplyr::select(x, !!f_enquo, !!coe_enquo)

  # Remove NAs ----------------
  df_NA <- is.na(df)
  cols_NA <- colSums(df_NA)
  rows_NA <- rowSums(df_NA)

  if (any(df_NA)){
    .msg_info("stat_lda_prepare: dropping NAs") # do we need to refine with where (?)
    df <- stats::na.omit(df)
  }

  # pulls f and drop useless levels
  f_naked <- dplyr::transmute_at(df, 1, forcats::as_factor)

  # check now for at least two
  if (nlevels(dplyr::pull(f_naked))<2)
    stop("stat_lda_prepare: f must have at least two levels")

  # numeric side now ----------
  dfn <- dplyr::select(df, -1)

  # handles constant
  cols_constant <- dfn %>%
    dplyr::summarise_all(stats::var) %>%
    dplyr::select_if(~.x < constant_below_var) %>% names()
  #.msg here
  if (length(cols_constant)>0){
    paste0("stat_lda_prepare: ",
           paste(cols_constant, collapse=", "),
           " were removed (var <", constant_below_var, ")") %>% .msg_info()
    dfn <- dfn %>% dplyr::select(-cols_constant)
  }

  # handles collinear ---------
  cormat <- dfn %>% stats::cor()
  cormat[upper.tri(cormat, diag = TRUE)] <- NA
  cols_collinear <- (cormat > collinear_above_cor) %>% apply(1, any, na.rm=TRUE) %>% which() %>% names()
  #.msg here
  if (length(cols_collinear)>0){
    paste0("stat_lda_prepare: ",
           paste(cols_collinear, collapse=", "),
           " were removed (cor > ", collinear_above_cor, ")") %>% .msg_info()
    dfn <- dfn %>% dplyr::select(-tidyselect::all_of(cols_collinear))
  }

  # we can affort a long one now
  coe_naked <- dfn

  # return these beauties
  list(df                = dplyr::bind_cols(f_naked, coe_naked),
       f_naked           = dplyr::pull(f_naked),
       coe_naked         = coe_naked,
       cols_constant     = cols_constant,
       cols_collinear    = cols_collinear,
       cols_NA           = cols_NA,
       rows_NA           = rows_NA)
}

# stat_lda0 -----------------------------------------------

#' Linear discriminant analysis
#'
#' Calculations are delegated to [MASS::lda]
#'
#' @param x `matrix or [tibble][tibble::tibble-package] containing the explanatory variables
#' forwarded to [MASS::lda]'s `x`
#' @param f `factor` forwarded to [MASS::lda]'s `grouping`. If missing takes first columns and remove it from x.
#' @param ... addtional parameters forwardes to [MASS::lda]
#' @param full `logical` whether to prepare useful components (default to `FALSE`)
#'
#' @details With `full=TRUE`, it is at least 6 times slower.
#'
#' @name stat_lda
#' @family  lda stats
#' @examples
#' (x <- dummy_df %>% stat_lda_prepare(foo2_NA, a:e))
#'
#' stat_lda0(x$coe_naked, x$f_naked)
#'
#' stat_lda(dummy_df, foo2_NA, a:e)
NULL

#' @describeIn stat_lda Vanilla lda
#' @export
stat_lda0 <- function(x, f, full=FALSE, ...){
  # delegate calculations to MASS
  mod <- MASS::lda(x, f, ...)
  mod_pred <- stats::predict(mod, x)

  # avoid a second call to lda
  # deduce the class ourselves
  f1 <- colnames(mod_pred$posterior)[apply(mod_pred$posterior, 1, which.max)] %>%
    factor(levels=levels(f))


  posterior_tbl <- tibble::tibble(actual    = f,
                                  predicted = f1,
                                  correct   = f==f1,
                                  posterior = apply(mod_pred$posterior, 1, max)) %>%
    dplyr::arrange(.data$actual, .data$predicted)

  # early return if one just want the prediction ==========
  if (!full)
    return(posterior_tbl)

  # otherwise return all details ==========================
  #
  # cross-validation section ------------------------------
  CV           <- posterior_tbl %>% dplyr::select(1:2) %>% table()

  # may exist a shorter path but this ensure all rows are returned
  CV_tbl <- CV %>% .CV_tbl

  # shortcuts vectors
  CV_accuracy <- sum(diag(CV))/sum(CV)
  CV_class_accuracy <- CV_tbl %>%
    dplyr::select(.data$actual, .data$prop_class) %>%
    dplyr::distinct() %>%
    tibble::deframe()

  # premultiplying matrix section -------------------------
  n <- nrow(x)
  lm.mod <- stats::lm(as.matrix(x) ~ f)
  dfw <- n - nlevels(f)
  SSw <- stats::var(lm.mod$residuals) * (n - 1)
  VCVw <- SSw/dfw
  LDs <- VCVw %*% mod$scaling

  # trace section -----------------------------------------
  trace_tbl <- tibble::tibble(axis           = colnames(mod_pred$x),
                              trace          = mod$svd^2,
                              proportion     = trace/sum(trace),
                              cum_proportion = cumsum(trace/sum(trace))) # onthe fly eg cumsum(proportion) R cmd check

  # return this beauty
  list(x                 = x,
       f                 = f,
       mod               = mod,
       scores            = tibble::as_tibble(mod_pred$x),
       posterior         = tibble::as_tibble(mod_pred$posterior),
       posterior_tbl     = posterior_tbl,
       CV                = CV,
       CV_tbl            = CV_tbl,
       CV_accuracy       = CV_accuracy,
       CV_class_accuracy = CV_class_accuracy,
       LDs = LDs) %>%
    structure(class=c("stat_lda", "list"))
}

# !full is 6 times faster

# microbenchmark::microbenchmark(stat_lda0(x, f),
# stat_lda0(x, f, TRUE), times=20)
# Unit: milliseconds
# expr                    min        lq        mean      median        uq
# stat_lda0(x, f)          3.243079  3.451461  3.828891  3.76649  4.195635
# stat_lda0(x, f, TRUE)    20.168418 21.801200 43.191688 22.52950 45.061394

# stat_lda ------------------------------------------------
#' @describeIn stat_lda Wrapped lda
#' @export
stat_lda <- function(x, f, ...){
  # tidyeval
  f_enquo   <- enquo(f)  ### todo: handles for formula here ??

  if (missing(...)) {

    if (Momocs2::coe_nb(x)==0)
      paste0("stat_lda: no coe column, and none column passed to '...'") %>% stop()

    coe <- Momocs2::coe_names(x)
    paste0("stat_lda: working on ", paste(coe, collapse=", ")) %>% .msg_info()
    coe_enquo <- rlang::expr(tidyselect::all_of(coe))
  } else {
    coe_enquo <- rlang::expr(c(...))
  }

  ready <- x %>% # add removed columns ?
    stat_lda_prepare(f=!!f_enquo, !!coe_enquo)

  res <- stat_lda0(x = ready$coe_naked, f=ready$f_naked, full=TRUE)
  res <- append(res, ready[c("posterior_tbl", "cols_constant", "cols_collinear")])
  class(res) <- c("stat_lda", "list")
  res
}

# print.stat_lda_full -------------------------------------
# x is typically a factor
.digest_balanced <- function(x){
  y <- table(x)
  if (length(unique(y))==1)
    return("balanced")
  else
    glue::glue_data(y, "unbalanced (N ranges from {min(y)} to {max(y)})")
}
# c("a", "a", "b", "b") %>% .digest_balanced()
# c("a", "a", "b") %>% .digest_balanced()

# x is a named vector
.digest_class_acc <- function(x){
  # might be shorter
  mini <- signif(x[x==min(x)], 3)
  medi <- signif(x[x==stats::median(x)], 3)
  maxi <- signif(x[x==max(x)], 3)
  glue::glue("min: {.digest_named_vec(mini)}, median: {.digest_named_vec(medi)}, max: {.digest_named_vec(maxi)}")
}

# prints nicely a named vector
# possibly with more than one name but a single value is expected
.digest_named_vec <- function(x){
  glue::glue("{as.numeric(x[1])} ({paste(names(x), collapse=', ')})")
}

# a direct alternative when names do not matter
.digest_numeric <- function(x){
  q <- stats::quantile(x, probs=c(0, 0.5, 1)) %>% as.numeric() %>% signif(3)
  glue::glue("min: {q[1]}, median: {q[2]}, max: {q[3]}")
}

#' @export
print.stat_lda <- function(x, ...){
  cli::cli_h1("Linear Discriminant Analysis")
  glue::glue("
    - {nrow(x$x)} observations
    - {ncol(x$x)} variables
    - {nlevels(x$f)} levels, {.digest_balanced(x$f)}") %>% cat()
  cli::cat_line()
  cli::cat_line()
  glue::glue("
    - Accuracy:            {signif(x$CV_accuracy, 3)}
    - Within classes:      {x$CV_class_accuracy %>% .digest_class_acc()}
    - Posterior (correct): {x$posterior_tbl %>% dplyr::filter(correct) %>% dplyr::pull(posterior) %>% .digest_numeric()}
                  ") %>% cat()
}


# stat_swing_lda --------


