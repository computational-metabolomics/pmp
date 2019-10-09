context("test-glog_transformation")

test_that("Glog function returns expected output", {
  out <- mv_imputation(df=testData$data, method="knn")
  out <- glog_transformation (df=out, classes=testData$class,
    qc_label="QC")
  expect_equal (out@scaled_peak_matrix, testData$glog_transformation)
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
      qc_label="QC")@scaled_peak_matrix[3,], regexp='Error!Lambda')
    testthat::expect_true(all(out[1:3] == out[4:6]))
})

test_that("glog function returns optimised lambda value if requested", {
  out <- mv_imputation(df=testData$data, method="knn")
  out <- glog_transformation (df=out, classes=testData$class,
    qc_label="QC")
  lambda <- as.integer(7865716)
  #testthat::expect_true(length(out) == 5)
  testthat::expect_true(lambda == as.integer(out@lambda_summary$lambda))
})

test_that("glog function returns plot of lambda optimisation if requested", {
  data <- mv_imputation(df=testData$data, method='knn')
  out <- glog_transformation (df=data, classes=testData$class, qc_label='QC')
  expect_equal(out@lambda_optimisation_plot[[1]][[9]]$x, "lambda")
  expect_equal(out@lambda_optimisation_plot[[1]][[9]]$y, "SSE")
})

