context("test-checkpeakmatrixorientation")

test_that("Output data matrix is in correct orientation", {
  out <- checkPeakMatrixOrientation(df=t(testData$data),classes = testData$class)
  expect_equal(out, testData$data)
})

test_that("Output data matrix is in correct orientation if no transpose is needed", {
  out <- checkPeakMatrixOrientation(df=testData$data,classes = testData$class)
  expect_equal(out, testData$data)
})

test_that("Function fails", {
  expect_error(checkPeakMatrixOrientation(df=testData$data, classes=testData$class[1:3]))
})

test_that("Function fails", {
  expect_error(checkPeakMatrixOrientation(df=testData$data[,1:6], classes=testData$class))
})