context ("test-filter_peaks_by_blank")

test_that ("Test that filter_peaks_by_blank returns expected output", {
  out <- filter_peaks_by_blank(df = testData$data, fold_change = 1.2, classes = testData$class,
                               blank_label = "Blank", qc_label = NULL, remove = FALSE, fraction_in_blank = 0)
  expect_equal(out, testData$filter_peaks_by_blank)
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
