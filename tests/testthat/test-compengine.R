# A unit test for compengine() function
if (require(testthat)) {
  context("Tests on input")
  test_that("tests for a non-vector object", {
    expect_that(compengine(matrix(0, 2, 2)), throws_error())
  })

  context("Tests on output")
  test_that("tests for compengine results on non-seasonal data", {
    z <- compengine(WWWusage)
    expect_equal(length(z), 16L)
    expect_equal(z[4], c(firstmin_ac = 21))
    expect_gt(z[5], 109.15)
    expect_gt(z[3], 0.27)
  })
  test_that("tests for compengine results on seasonal data", {
    z <- compengine(USAccDeaths)
    expect_that(length(z), equals(16L))
    expect_equal(z[4], c(firstmin_ac = 6))
    expect_gt(z[6], 1.83)
    expect_lt(z[3], -0.0647)
  })
  test_that("tests for compengine results on data with missing values", {
    y_WWWusage <- WWWusage
    y_WWWusage[c(16:17, 78)] <- NA
    z <- compengine(y_WWWusage)
    expect_equal(length(which(is.na(z))), 0)
    expect_gt(z[3], 0.2845)
    expect_equal(z[4], c(firstmin_ac = 21))
  })
}
