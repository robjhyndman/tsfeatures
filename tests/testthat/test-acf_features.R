# A unit tests for acf_features() function

if (require(testthat)) {
  context("Tests on output")
  test_that("test for acf_features() results on non-seasonal ts data", {
    z <- acf_features(WWWusage)
    expect_equal(length(z), 6L)
    expect_gt(z[1], 0.96)
    expect_gt(z[2], 4.19)
    expect_gt(z[3], 0.79)
    expect_gt(z[4], 1.40)
    expect_gt(z[5], 0.17)
    expect_gt(z[6], 0.33)
  })
  test_that("test for acf_features() results on seasonal ts data", {
    z <- acf_features(USAccDeaths)
    expect_equal(length(z), 7L)
    expect_gt(z[1], 0.70)
    expect_gt(z[2], 1.20)
    expect_gt(z[3], 0.02)
    expect_gt(z[4], 0.27)
    expect_gt(z[5], -0.49)
    expect_gt(z[6], 0.74)
    expect_gt(z[7], 0.62)
  })
}
