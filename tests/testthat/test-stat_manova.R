test_that(".cross works", {
  k <- 5
  x <- .cross(c(letters[1:k]))
  expect_is(x, "tbl")
  expect_equal(nrow(x), (5*4 / 2))
})
