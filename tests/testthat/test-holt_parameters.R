# A unit tests for holt_parameters() function

if (require(testthat)) {
  context("Tests on output")
  test_that("test for holt_parameters() results on non-seasonal ts data", {
    z <- holt_parameters(WWWusage)
    expect_equal(length(z), 2L)
    expect_gt(z[1], 0.99)
    expect_gt(z[2], 0.99)
  })
  test_that("test for holt_parameters() results on seasonal ts data", {
    z <- holt_parameters(USAccDeaths)
    expect_equal(length(z), 2L)
    expect_gt(z[1], 0.96)
    expect_gt(z[2], 0.00)
  })
}
