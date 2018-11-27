context("test-glog_transformation")

test_that("Glog function returns expected output", {
  out <- mv_imputation(df=testData$data, method = "knn")
  expect_equal (glog_transformation (df=out, classes = testData$class, qc_label = "QC"), testData$glog_transformation)
})

test_that("Glog function fails if qc_label is wrong or QC samples don't exist", {
  out <- mv_imputation(df=testData$data, method = "knn")
  expect_error (glog_transformation (df=out, classes = testData$class, qc_label = "A"))
})
