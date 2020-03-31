test_that(".CV_tbl works", {
  set.seed(2329)
  m <- matrix(sample(1:100), 10)
  cv <- .CV_tbl(m)
  expect_is(cv, "tbl")
  expect_equal(nrow(cv), prod(dim(m)))
  expect_equal(colnames(cv), c("actual", "predicted", "n", "prop", "n_class", "prop_class", "n_total"))
})


test_that("stat_lda_prepare works", {
  # a nice tibble with many problems

  set.seed(2329)
  k=100

  df <- tibble::tibble(
    # q variables, some collinear,
    a=rnorm(k), b=jitter(a), c=runif(k),
    # now constants, collinear or not
    d=pi, e=jitter(d), f=3,
    # factors now
    foo1=factor("a"),
    foo2=sample(LETTERS[1:12], size=k, replace = TRUE) %>% factor(),
    foo3=sample(c("yes", "no"), size=k, replace=TRUE) %>% factor(),
  )
  e_NA <- df$e
  e_NA[c(5, 12)] <- NA

  foo2_NA <- df$foo2
  foo2_NA[c(7, 34)] <- NA

  df <- df %>% dplyr::mutate(e_NA=e_NA, foo2_NA=foo2_NA)

  expect_message(df_prepared <- df %>% stat_lda_prepare(foo2_NA, a:f))
  expect_is(df_prepared, "list")

  # no NAs
  expect_equal(sum(is.na(df_prepared$df)), 0)

  # no more constant
  variances <- df_prepared$df %>% dplyr::select(-1) %>% dplyr::summarise_all(stats::var) %>% unlist()
  expect_true(all(variances > 1e-5))

  # no more collinear
  correlations <- df_prepared$df %>% dplyr::select(-1) %>% stats::cor()
  correlations[upper.tri(correlations, diag=TRUE)] <- 0
  expect_equal(sum(correlations>(1-1e-5)), 0)

  # now it should be ok so no messaging
  expect_silent(df_prepared2 <- df_prepared$df %>% stat_lda_prepare(foo2_NA, a, c))
  expect_equal(df_prepared$df %>% dplyr::select(1, a, c), df_prepared2$df)
})

