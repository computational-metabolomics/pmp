context("test data input/output functions")

test_that("Output data matrix is in correct orientation", {
  out <- check_peak_matrix(df=t(testData$data), classes=testData$class)
  expect_equal(out, testData$data)
})

test_that("Output data matrix is in correct orientation if no transpose is needed", {
  out <- check_peak_matrix(df=testData$data, classes=testData$class)
  expect_equal(out, testData$data)
})

test_that("Function fails", {
  expect_error(pmp:::check_peak_matrix(df=testData$data, classes=testData$class[1:3]))
})

test_that("Function fails", {
  expect_error(check_peak_matrix(df=testData$data[,1:6], classes=testData$class))
})

test_that("Function fails, input is character", {
  out <- testData$data
  out[ ,1] <- as.character(out[ ,1])
  expect_error(pmp:::check_peak_matrix(df=out, classes=testData$class))
})

test_that("Function works if class labels are not provided", {
  out <- check_peak_matrix(df=testData$data)
  expect_equal (out, testData$data)
})

test_that("Function works if class labels are not provided and matrix needs to be transposed", {
  expect_warning (out <- check_peak_matrix(df=t(testData$data)))
  expect_equal(out, testData$data)
})

test_that("Function works if class labels are not provided and matrix needs to be transposed,
          and data type is retained", {
  df <- data.frame(t(testData$data))
  expect_warning (out <- check_peak_matrix(df=df))
})

test_that("Function returns warning if input peak matrix has the same number of samples and features", {
  expect_warning (out <- check_peak_matrix(df=testData$data[1:9, ], 
    classes=testData$class))
})

test_that("return_original_data_structure work with bioconductor DataFrame object", {
  df <- S4Vectors::DataFrame(A=c(1:5), B=c(1:5))
  df <- pmp:::check_input_data (df)
  meta_data <- metadata(df)
  meta_data$processing_history <- "testthat"
  metadata(df) <- meta_data
  expect_silent(df <- pmp:::return_original_data_structure(df))
})

test_that("SummarizedExperiment has original data structure flag added", {
  df <- MTBLS79[, MTBLS79$Batch == 1]
  expect_silent(df <- check_input_data(df))
  expect_true(metadata(df)$original_data_structure == "SummarizedExperiment")
  expect_silent(df <- pmp:::return_original_data_structure(df))
  expect_true (class(df)[1] == "SummarizedExperiment")
})
