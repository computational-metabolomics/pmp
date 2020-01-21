context ("test-filter_peaks_by_blank")

test_that ("Test that filter_peaks_by_blank returns expected output", {
  out <- filter_peaks_by_blank(df=testData$data, fold_change=1.2, 
    classes=testData$class, blank_label="Blank", qc_label=NULL, 
    remove_samples=FALSE, remove_peaks=TRUE, fraction_in_blank=0)
  attributes(out)$processing_history <- NULL
  expect_equal(as.matrix(attributes(out)$flags),
    as.matrix(testData$filter_peaks_by_blank$flags))
  attributes(out)$flags <- NULL
  expect_equal(as.matrix(out), testData$filter_peaks_by_blank$df)
  
})

test_that ("Test that filter_peaks_by_blank returns expected output with QC samples specified", {
  out <- filter_peaks_by_blank(df=testData$data, fold_change=1.2, 
    classes=testData$class, blank_label="Blank", qc_label="QC",
    remove_samples=FALSE, remove_peaks=TRUE, fraction_in_blank=0)
  expect_equal(processing_history(out)[[1]][[1]], 1.2)
    attributes(out)$processing_history <- NULL
    attributes(out)$flags <- NULL
  expect_equal(out, testData$filter_peaks_by_blank_qc$df)
})

test_that ("Test that filter_peaks_by_blank removes blank samples, but keeps features", {
  out <- filter_peaks_by_blank(df=testData$data, fold_change=1.2, 
    classes=testData$class, blank_label="Blank", qc_label=NULL, 
    remove_samples=TRUE, remove_peaks=FALSE,fraction_in_blank=0)
  expect_equal(as.matrix(attributes(out)$flags), 
    as.matrix(testData$filter_peaks_by_blank_remove_blanks$flags))
    attributes(out)$flags <- NULL
  expect_equal(out, testData$filter_peaks_by_blank_remove_blanks$df)
  
})

context ("test-filter_peaks_by_fraction")

test_that ("Test that filter_peaks_by_fraction returns expected output, method=QC", {
  out <- filter_peaks_by_fraction(df=testData$data, min_frac=1, 
  classes=testData$class, method="QC", qc_label="QC")
  attributes(out)$processing_history <- NULL
  expect_equal(attributes(out)$flags, testData$filter_peaks_by_fraction$flags)
  attributes(out)$flags <- NULL
  expect_equal(out, testData$filter_peaks_by_fraction$df)
  
})

test_that ("Test that filter_peaks_by_fraction returns expected output, method=within", {
  out <- filter_peaks_by_fraction(df=testData$data, min_frac=1, 
  classes=testData$class, method="within", qc_label="QC")
  attributes(out)$processing_history <- NULL
  expect_equal(attributes(out)$flags, testData$filter_peaks_by_fraction_within$flags)
  attributes(out)$flags <- NULL
  expect_equal(out, testData$filter_peaks_by_fraction_within$df)
  
})

test_that ("Test that filter_peaks_by_fraction returns expected output, method=across", {
  out <- filter_peaks_by_fraction(df=t(testData$data), min_frac=0.6, 
    classes=testData$class, method="across", qc_label="QC")
  attributes(out)$processing_history <- NULL
  expect_equal(attributes(out)$flags, testData$filter_peaks_by_fraction_across$flags)
  attributes(out)$flags <- NULL
  expect_equal(out, testData$filter_peaks_by_fraction_across$df)
})

context ("test-filter_peaks_by_rsd")

test_that("Test that filter_peaks_by_rsd returns expected output", {
  out <- filter_peaks_by_rsd(df=t(testData$data), max_rsd=20, 
    classes=testData$class, qc_label="QC")
  attributes(out)$processing_history <- NULL
  expect_equal (attributes(out)$flags, testData$filter_peaks_by_rsd$flags)
  attributes(out)$flags <- NULL
  expect_equal (out, testData$filter_peaks_by_rsd$df)
  
})

context ("test-filter_samples_by_mv")

test_that("Test that filter_samples_by_mv returns expected output", {
  
  out <- testData$data
  out[c(1:25),c(2,9)] <- NA
  expect_warning(out <- filter_samples_by_mv (df=t(out), max_perc_mv=0.8))
  attributes(out)$processing_history <- NULL
  flags <- attributes(out)$flags
  colnames(flags) <- c("perc_mv", "flags")
  expect_equal (flags, testData$filter_samples_by_mv$flags)
  attributes(out)$flags <- NULL
  expect_equal (out, testData$filter_samples_by_mv$df)
  
})

test_that("Remove peaks filter returns expecteed output", {
  
  rem_index <-testData$remove_peaks$rem_index
  out <- remove_peaks(df=testData$data, rem_index=rem_index)
  attributes(out)$processing_history <- NULL
  expect_equal (out, testData$remove_peaks$rem_output)
  
})

test_that("Remove peaks filter fails if input is not logical vector", {
  
  rem_index <- as.factor(testData$remove_peaks$rem_index)
  expect_error(out <- remove_peaks(df=testData$data, rem_index=rem_index))
  
})

test_that("Blank filter fails if one column is a factor", {
  data <- as.data.frame(testData$data)
  data$QC_2 <- as.factor(data$QC_2)
  testthat::expect_error(filter_peaks_by_blank(df=data, 
    fold_change=1.2, classes=testData$class, 
    blank_label="Blank", qc_label=NULL, remove=FALSE, 
    fraction_in_blank=0))
})

test_that("Blank filter fails if all data are data type 'character'", {
  data <- as.character(testData$data)
  testthat::expect_error(filter_peaks_by_blank(df=data, 
    fold_change=1.2, classes=testData$class, blank_label="Blank", 
    qc_label=NULL, remove=FALSE, fraction_in_blank=0))
})
