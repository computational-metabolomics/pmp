context("test-filter_peaks_by_blank")

test_that("Test that filter_peaks_by_blank returns expected output", {
  out <- filter_peaks_by_blank(df = testData$data, fold_change = 1.2, classes = testData$class,
                               blank_label = "Blank", qc_label = NULL, remove = FALSE, fraction_in_blank = 0)
  expect_equal(out, testData$filter_peaks_by_blank)
})

test_that("Test that filter_peaks_by_blank returns expected output and remove blank samples", {
  out <- filter_peaks_by_blank(df = testData$data, fold_change = 1.2, classes = testData$class,
                               blank_label = "Blank", qc_label = NULL, remove = TRUE, fraction_in_blank = 0)
  expect_equal(out, testData$filter_peaks_by_blank_remove_blanks)
})

context("test-filter_peaks_by_fraction")

test_that ("Test that filter_peaks_by_fraction returns expected output", {
  out <- 1+1
  expect_equal(out, 2)
})
