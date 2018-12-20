# A unit tests for pacf_features() function

if (require(testthat)) {
  context("Tests on output")
  test_that("test for pacf_features() results on non-seasonal ts data", {
    z <- pacf_features(WWWusage)
    expect_equal(length(z), 3L)
    expect_gt(z[1], 1.03)
    expect_gt(z[2], 0.80)
    expect_gt(z[3], 0.22)
  })
  test_that("test for pacf_features() results on seasonal ts data", {
    z <- pacf_features(USAccDeaths)
    expect_equal(length(z), 4L)
    expect_gt(z[1], 0.63)
    expect_gt(z[2], 0.09)
    expect_gt(z[3], 0.38)
    expect_gt(z[4], 0.12)
  })
}
