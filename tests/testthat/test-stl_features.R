# A unit test for stl_features() function
if (require(testthat)) {
  context("Tests on input")
  test_that("tests for a non-vector object", {
    expect_that(stl_features(matrix(0, 2, 2)), throws_error())
  })

  context("Tests on output")
  test_that("tests for stl_feature results on non-seasonal data", {
    z <- stl_features(WWWusage)
    expect_equal(length(z), 8L)
    expect_equal(z[1], c(nperiods = 0))
    expect_equal(z[2], c(seasonal_period = 1))
    expect_gt(z[3], 0.98)
  })
  test_that("tests for stl_feature results on seasonal ts data", {
    z <- stl_features(USAccDeaths)
    expect_that(length(z), equals(11L))
    expect_equal(z[1], c(nperiods = 1))
    expect_equal(z[2], c(seasonal_period = 12))
    expect_gt(z[3], 0.78)
  })
  test_that("tests for stl_feature results on seasonal msts data", {
    z <- stl_features(forecast::taylor)
    expect_that(length(z), equals(15L))
    expect_equal(z[1], c(nperiods = 2))
    expect_equal(z[2], c(seasonal_period1 = 48))
    expect_equal(z[3], c(seasonal_period2 = 336))
    expect_gt(z[4], 0.79)
  })
}
