context("test-glog_transformation")

test_that("Glog function returns expected output", {
  out <- mv_imputation(df=testData$data, method="knn")
  out <- glog_transformation (df=out, classes=testData$class,
    qc_label="QC")
  expect_equal (assay(out), as.matrix(testData$glog_transformation))
  
  lambda <- as.integer(7865716)
  testthat::expect_true(lambda == as.integer(metadata(out)$glog_scaling$lambda))
  
  #out <- getGlogLambdaOptimisationPlot(out)
  
  #expect_equal(out[[9]]$x, "lambda")
  #expect_equal(out[[9]]$y, "SSE")
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
    out <- assay(out)[3,]
    testthat::expect_true(all(out[1:3] == out[4:6]))
})
