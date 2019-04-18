context("test-normalise_to_total_sum")

test_that("normalise to sum returns correct output", {
  expect_equal(normalise_to_sum (testData$data), testData$normalise_to_sum)
})

test_that("normalise to sum returns correct output when matrix needs to be transposed", {
  expect_warning (out <- normalise_to_sum (t(testData$data)))
  expect_equal(out, testData$normalise_to_sum)
})

context("test-pqn_normalisation")

test_that("PQN normalisation returns correct output", {
  out <- pqn_normalisation(df=testData$data, classes=testData$class, qc_label = "QC")
  expect_equal(out, testData$pqn_normalisation)
})

test_that("PQN normalisation returns correct output, when all samples are used to calculate correction factor", {
  out <- pqn_normalisation(df=testData$data, classes=testData$class, qc_label = "all")
  expect_equal(out, testData$pqn_normalisation_all)
})

test_that("PQN normalisation returns correct output when matrix needs to be transposed", {
  out <- pqn_normalisation(df=t(testData$data), classes=testData$class, qc_label = "QC")
  expect_equal(out, testData$pqn_normalisation)
})