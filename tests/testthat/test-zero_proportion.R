# A unit tests for zero_proportion() function

if (require(testthat)) {
  context("Tests on output")
  test_that("test for zero_proportion() results on ts data", {
    z <- as.ts(c(12, 13, 4, 0, 0.0, 23, 45, 67, 89))
    y <- expect_equal(length(z), 1L)
    expect_gt(y, 2)
  })

}
