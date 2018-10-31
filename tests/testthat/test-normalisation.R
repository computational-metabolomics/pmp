context("test-normalise_to_total_sum")

test_that("normalise to sum returns correct output", {
  expect_equal(normalise_to_sum (testData$data), testData$normalise_to_sum)
})

test_that("normalise to sum returns correct output when matrix needs transposing", {
  expect_warning (out <- normalise_to_sum (t(testData$data)))
  expect_equal(out, testData$normalise_to_sum)
})