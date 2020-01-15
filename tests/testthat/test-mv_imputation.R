context("test-mv_imputation")

test_that("Missing value imputation using knn method returns expected output", {
  out <- mv_imputation(df=testData$data, method="knn", check_df=FALSE)
  attributes(out)$processing_history <- NULL
  out <- as.data.frame(out)
  expect_equal(out, testData$mv_imputation_knn)
})

test_that("Missing value imputation using rf method returns expected output", {
  expect_warning(out <- mv_imputation(df=testData$data, method="rf"))
  attributes(out)$processing_history <- NULL
  out <- as.data.frame(out)
  expect_equal(out, testData$mv_imputation_rf)
})

test_that("Missing value imputation using bpca method returns expected output", {
  expect_warning (out <- mv_imputation(df=t(testData$data), method="bpca"))
  attributes(out)$processing_history <- NULL
  out <- as.data.frame(out)
  expect_equal(out, testData$mv_imputation_bpca)
})

test_that("Missing value imputation using sv method returns expected output", {
  expect_warning (out <- mv_imputation(df=t(testData$data), method="sv"))
  attributes(out)$processing_history <- NULL
  out <- as.data.frame(out)
  expect_equal(out, testData$mv_imputation_sv)
})

test_that("Missing value imputation using mn method returns expected output", {
  expect_warning (out <- mv_imputation(df=t(testData$data), method="mn"))
  attributes(out)$processing_history <- NULL
  out <- as.data.frame(out)
  expect_equal(out, testData$mv_imputation_mn)
})

test_that("Missing value imputation using md method returns expected output", {
  expect_warning (out <- mv_imputation(df=t(testData$data), method="md"))
  attributes(out)$processing_history <- NULL
  out <- as.data.frame(out)
  expect_equal(out, testData$mv_imputation_md)
})

test_that("Missing value imputation stops if all values in the row are NA's", {
    out <- testData$data
    out[3,] <- NA
    expect_error(out <- mv_imputation(df=out, method="md"))
})

test_that("Missing value imputation stops if all values in the column are NA's", {
  out <- testData$data
  out[,3] <- NA
  expect_error(out <- mv_imputation(df=out, method="md"))
})

test_that("Missing value imputation stops if wrong method is selected", {
  out <- testData$data
  expect_error(out <- mv_imputation(df=out, method="Kn"))
})
