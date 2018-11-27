context("test-dopca")

test_that("dopca creates correct plot object", {
  out <- doPCA (Data = testData$prepareData_mv_impute_glog_pqn$Data, classes = testData$prepareData_mv_impute_glog_pqn$classes, 
                qc_label = "QC", PQN=T, glogScaling=T, scale=F, mv_impute = T, loadings = T)
  expect_equal (out, testData$doPCA)
})

test_that("dopca creates correct plot object without sample labels", {
  out <- doPCA (Data = testData$prepareData_mv_impute_glog_pqn$Data, classes = testData$prepareData_mv_impute_glog_pqn$classes, 
                qc_label = "QC", PQN=T, glogScaling=T, scale=F, mv_impute = T, loadings = F, labels = "none")
  expect_equal (out, testData$doPCA_no_labels)
})