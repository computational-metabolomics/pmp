context("test-glog_transformation")

test_that("Glog function returns expected output", {
  out <- mv_imputation(df=testData$data, method="knn")
  out <- glog_transformation (df=out, classes=testData$class,
    qc_label="QC")
  
  lambda <- as.integer(7865716)
  testthat::expect_true(lambda == as.integer(attributes(out)$
    processing_history$glog_transformation$lambda_opt))
  
  expect_equal (as.data.frame(out), testData$glog_transformation)
})

test_that("Glog function fails if qc_label is wrong or QC samples don't exist", {
  out <- mv_imputation(df=testData$data, method="knn")
  expect_error (glog_transformation (df=out, classes=testData$class,
    qc_label="A"))
})

test_that("Scale factor is applied if optimisation fails", {
    df <- iris
    classes <- iris$Species
    df <- df[,-5]
    expect_output(glog_transformation(df,classes,qc_label='versicolor'),
     regexp='Error!Lambda')
})

test_that("If feature variance is 0 replace it with small value", {
    out <- mv_imputation(df=testData$data, method="knn")
    out[3,] <- 3
    expect_output(out <- glog_transformation (df=out, classes=testData$class,
      qc_label="QC"), regexp='Error!Lambda')
    out <- out[3,]
    testthat::expect_true(all(out[1:3] == out[4:6]))
})
