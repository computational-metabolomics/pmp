context("test-glog_transformation")

test_that("Glog function returns expected output", {
  out <- mv_imputation(df=t(testData$data), method = "knn")
  expect_equal (glog_transformation (df=out, classes = testData$class, qc_label = "QC"), testData$glog_transformation)
})
