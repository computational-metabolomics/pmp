context("test-glog_transformation")

test_that("Glog function returns expected output", {
  out <- mv_imputation(df=testData$data, method="knn")
  expect_equal (glog_transformation (df=out, classes=testData$class,
    qc_label="QC"), testData$glog_transformation)
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
      qc_label="QC")[3,], regexp='Error!Lambda')
    testthat::expect_true(all(out[1:3] == out[4:6]))
})

test_that("glog function returns optimised lambda value if requested", {
  out <- mv_imputation(df=testData$data, method="knn")
  out <- glog_transformation (df=out, classes=testData$class,
    qc_label="QC", store_lambda=TRUE)
  lambda <- as.integer(7865716)
  testthat::expect_true(length(out) == 4)
  testthat::expect_true(lambda == as.integer(out[[2]]))
})

test_that("glog function returns plot of lambda optimisation if requested", {
  data <- mv_imputation(df=testData$data, method='knn')
  out <- glog_transformation (df=data, classes=testData$class, qc_label='QC',
                              store_lambda=TRUE)
  data_qc <- data[,testData$class=="QC"]
  expect_silent(g <- glog_plot_optimised_labmda (optimised_lambda = out[[2]], data_qc = data_qc,
      upper_lim=10^9))
  expect_equal(g[[9]]$x, "lambda")
  expect_equal(g[[9]]$y, "SSE")
})

