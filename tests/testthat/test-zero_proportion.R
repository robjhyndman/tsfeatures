# A unit tests for zero_proportion() function

if (require(testthat)) {
  context("Tests on output")
  test_that("test for zero_proportion() ", {
    z <- zero_proportion(as.ts(c(0, 0, 3, 1, 2, 0)))
    expect_equal(length(z), 1L)
    expect_equal(z[1], 0.5)
  })
}
