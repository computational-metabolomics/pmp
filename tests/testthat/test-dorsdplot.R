context("test-dorsdplot")

test_that("That violin plot of RSD values is created", {
  out <- doRSDplot (RSD = testData$doRSD, plotTitle = "Test", subtitle = "sub title")
  expect_equal(out, testData$doRSDplot)
})
