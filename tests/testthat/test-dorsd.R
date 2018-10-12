context("test-dorsd")

test_that("Check sonsistency of doRSD function output", {
  out <- doRSD (testData$data, testData$class)  
  expect_equal(out , testData$doRSD)
})
