context ("test-filter_peaks_by_blank")

test_that ("Test that filter_peaks_by_blank returns expected output", {
  out <- filter_peaks_by_blank(df = testData$data, fold_change = 1.2, classes = testData$class,
                               blank_label = "Blank", qc_label = NULL, remove = FALSE, fraction_in_blank = 0)
  expect_equal(out, testData$filter_peaks_by_blank)
})

test_that ("Test that filter_peaks_by_blank returns expected output with QC samples specified", {
  out <- filter_peaks_by_blank(df = testData$data, fold_change = 1.2, classes = testData$class,
                               blank_label = "Blank", qc_label = "QC", remove = FALSE, fraction_in_blank = 0)
  expect_equal(out, testData$filter_peaks_by_blank_qc)
})

test_that ("Test that filter_peaks_by_blank returns expected output and remove blank samples", {
  out <- filter_peaks_by_blank(df = testData$data, fold_change = 1.2, classes = testData$class,
                               blank_label = "Blank", qc_label = NULL, remove = TRUE, fraction_in_blank = 0)
  expect_equal(out, testData$filter_peaks_by_blank_remove_blanks)
})

context ("test-filter_peaks_by_fraction")

test_that ("Test that filter_peaks_by_fraction returns expected output, method=QC", {
  out <- filter_peaks_by_fraction(df = testData$data, min_frac = 1, classes = testData$class, method = "QC", qc_label = "QC")
  expect_equal(out, testData$filter_peaks_by_fraction)
})

test_that ("Test that filter_peaks_by_fraction returns expected output, method=within", {
  out <- filter_peaks_by_fraction(df = testData$data, min_frac = 1, classes = testData$class, method = "within", qc_label = "QC")
  expect_equal(out, testData$filter_peaks_by_fraction_within)
})

test_that ("Test that filter_peaks_by_fraction returns expected output, method=across", {
  out <- filter_peaks_by_fraction(df = t(testData$data), min_frac = 0.6, classes = testData$class, method = "across", qc_label = "QC")
  expect_equal(out, testData$filter_peaks_by_fraction_across)
})

context ("test-filter_peaks_by_rsd")

test_that("Test that filter_peaks_by_rsd returns expected outpu", {
  out <- filter_peaks_by_rsd(df=t(testData$data), max_rsd = 20, classes = testData$class, qc_label = "QC")
  expect_equal (out, testData$filter_peaks_by_rsd)
})

context ("test-filter_peaks_by_rsd")

test_that("Test that filter_peaks_by_rsd returns expected outpu", {
  out <- filter_peaks_by_rsd(df=t(testData$data), max_rsd = 20, classes = testData$class, qc_label = "QC")
  expect_equal (out, testData$filter_peaks_by_rsd)
})

context ("test-filter_samples_by_mv")

test_that("Test that filter_peaks_by_rsd returns expected output", {
  
  out <- testData$data
  out[c(1:25),c(2,9)] <- NA
  expect_warning(out <- filter_samples_by_mv (df = t(out), max_perc_mv = 0.8))
  expect_equal (out, testData$filter_samples_by_mv)
})

test_that("Remove peaks filter returns expecteed output", {
  
  rem_index <- testData$remove_peaks$rem_index
  out <- remove_peaks(df=testData$data, rem_index = rem_index)
  expect_equal (out, testData$remove_peaks$rem_output)
  
})

test_that("Remove peaks filter fails if input is not logical vector", {
  
  rem_index <- as.factor(testData$remove_peaks$rem_index)
  expect_error(out <- remove_peaks(df=testData$data, rem_index = rem_index))
  
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
                                               fold_change=1.2, classes=testData$class, 
                                               blank_label="Blank", qc_label=NULL, remove=FALSE, 
                                               fraction_in_blank=0))
})

