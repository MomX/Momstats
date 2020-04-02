test_that(".CV_tbl works", {
  set.seed(2329)
  m <- matrix(sample(1:100), 10)
  cv <- .CV_tbl(m)
  expect_is(cv, "tbl")
  expect_equal(nrow(cv), prod(dim(m)))
  expect_equal(colnames(cv), c("actual", "predicted", "n", "prop", "n_class", "prop_class", "n_total"))
})

# a nice tibble with many problems
df <- dummy_df

test_that("stat_lda_prepare works", {

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

test_that("stat_lda0 works fine", {
  df0 <- df %>% stat_lda_prepare(foo2_NA, a:f)
  # positionally
  lda1 <- stat_lda0(df0$coe_naked, df0$f_naked)
  expect_is(lda1, "tbl")
})

test_that("stat_lda works fine", {
  df <- dummy_df
  z <- df %>% stat_lda(foo2_NA, a:f)
  expect_is(z, "stat_lda")
  expect_output(print(z))
  # positionally
  expect_is(stat_lda(df, foo2_NA, a), "stat_lda")
  # error since foo1 single level
  expect_error(stat_lda(df, foo1, a:b))
  # this one should work, also add one more column
  expect_is(stat_lda(df, foo2, a:b, e_NA), "stat_lda")
})

test_that(".digest_* work", {
  expect_equal(c("a", "a", "b", "b") %>% .digest_balanced(), "balanced")
  expect_equal(c("a", "a", "b") %>% .digest_balanced(), "unbalanced (N ranges from 1 to 2)")

  expect_equal(1:10 %>% .digest_numeric(), "min: 1, median: 5.5, max: 10")
  expect_equal(1 %>% .digest_numeric(), "min: 1, median: 1, max: 1")

  class_acc <- c(plop=0.5, plip=1.2, yup=5)

  expect_equal(class_acc %>% .digest_class_acc(), "min: 0.5 (plop), median: 1.2 (plip), max: 5 (yup)")
  expect_equal(class_acc %>% .digest_named_vec(), "0.5 (plop, plip, yup)")
})


