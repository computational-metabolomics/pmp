context("test-preparedata")

test_that("Prepare data returns expected output, mv imputation only", {
  out <- prepareData(Data=testData$data, classes=testData$class,
              blank = "Blank", PQN=F, mv_impute = T, glogScaling = F,
              qc_label = "QC", ignorelabel = "Removed")
  expect_equal(out, testData$prepareData_mv_impute)
})

test_that("Prepare data returns expected output, mv imputation and glog", {
  out <- prepareData(Data=testData$data, classes=testData$class,
                     blank = "Blank", PQN=F, mv_impute = T, glogScaling = T,
                     qc_label = "QC", ignorelabel = "Removed")
  expect_equal(out, testData$prepareData_mv_impute_glog)
})

test_that("Prepare data returns expected output, mv imputation, glog and pqn", {
  out <- prepareData(Data=testData$data, classes=testData$class,
                     blank = "Blank", PQN=T, mv_impute = T, glogScaling = T,
                     qc_label = "QC", ignorelabel = "Removed")
  expect_equal(out, testData$prepareData_mv_impute_glog_pqn)
})

test_that("If all values are missing the row or columns will be removed", {
  Data <- testData$data
  Data[7, ] <- NA
  out <- prepareData(Data=Data, classes=testData$class,
                     blank = "Blank", PQN=F, mv_impute = T, glogScaling = F,
                     qc_label = NULL, ignorelabel = "Removed")
  expect_equal(out$Data, testData$prepareData_test_NA)
})

