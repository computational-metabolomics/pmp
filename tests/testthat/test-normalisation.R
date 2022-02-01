context("test-normalise_to_total_sum")

test_that("normalise to sum returns correct output", {
  out <- normalise_to_sum (testData$data)
  attributes(out)$processing_history <- NULL
  expect_equal(out, testData$normalise_to_sum)
})

test_that("normalise to sum returns correct output when matrix needs to be transposed", {
  expect_warning (out <- normalise_to_sum(t(testData$data)))
  attributes(out)$processing_history <- NULL
  expect_equal(out, testData$normalise_to_sum)
})

test_that("normalise to sum returns correct output when there are less features than samples", {
  out <- normalise_to_sum(testData$data[1:8,], check_df=FALSE)
  expect_true(all(round(apply (out, 2, sum, na.rm=TRUE), 0) == 100L))
})

context("test-pqn_normalisation")

test_that("PQN normalisation returns correct output", {
  out <- pqn_normalisation(df=testData$data, classes=testData$class, 
    qc_label="QC")
  attributes(out)$processing_history <- NULL
  expect_equal(as.vector(attributes(out)$flags),
    testData$pqn_normalisation$coef)
  attributes(out)$flags <- NULL
  expect_equal(out, testData$pqn_normalisation$df)
})

test_that("PQN normalisation returns correct output, when all samples are used to calculate correction factor", {
  out <- pqn_normalisation(df=testData$data, classes=testData$class, 
    qc_label="all")
  attributes(out)$processing_history <- NULL
  expect_equal(as.vector(attributes(out)$flags),
               testData$pqn_normalisation_all$coef)
  attributes(out)$flags <- NULL
  expect_equal(out, testData$pqn_normalisation_all$df)
})

test_that("PQN normalisation returns correct output when matrix needs to be transposed", {
  out <- pqn_normalisation(df=t(testData$data), classes=testData$class, 
    qc_label="QC")
  attributes(out)$processing_history <- NULL
  expect_equal(as.vector(attributes(out)$flags),
               testData$pqn_normalisation$coef)
  attributes(out)$flags <- NULL
  expect_equal(out, testData$pqn_normalisation$df)
})

test_that("PQN normalisation doesn't crash if only one feature is selected for normalisation", {
  out <- matrix(nrow=3, ncol=9)
  out <- rbind(out, testData$data[1, ])
  expect_warning(out <- pqn_normalisation(df=out, classes=testData$class, qc_label="QC"))
  expect_true(nrow(out) == 4)
})


test_that("PQN computation of reference works as expected for mean and median",{
  out <- pqn_normalisation(df=testData$data, classes=testData$class, 
    qc_label="QC",ref_method = 'mean')
  mean_ref=attributes(out)$processing_history$pqn_normalisation$computed_ref
  expect_equal(mean_ref[[1]],106308.238,tolerance = 0.0005)
  
  out <- pqn_normalisation(df=testData$data, classes=testData$class, 
    qc_label="QC",ref_method = 'median')
  median_ref=attributes(out)$processing_history$pqn_normalisation$computed_ref
  expect_equal(median_ref[[1]],88597.333,tolerance = 0.0005)
})


test_that("PQN check reference when appling additional filtering.",{
  out <- pqn_normalisation(df=testData$data, classes=testData$class, 
    qc_label="QC",qc_frac = 1,sample_frac=1,ref_method = 'mean')
  mean_ref=attributes(out)$processing_history$pqn_normalisation$computed_ref
  expect_equal(mean_ref[[1]],19855.36,tolerance = 0.0005)
  
  out <- pqn_normalisation(df=testData$data, classes=testData$class, 
    qc_label="QC",qc_frac = 1,sample_frac=1,ref_method = 'median')
  median_ref=attributes(out)$processing_history$pqn_normalisation$computed_ref
  expect_equal(median_ref[[1]],10773.31,tolerance = 0.0005)
})



test_that("PQN throws an error if incorrect ref method provided.",{
    expect_error(out <- pqn_normalisation(df=testData$data, classes=testData$class, 
        qc_label="QC",qc_frac = 1,sample_frac=1,ref_method = 'cake'))
})

test_that("PQN throws an error if QC filtering is too strict.",{
    expect_error(out <- pqn_normalisation(df=testData$data, classes=testData$class, 
        qc_label="QC",qc_frac = 1.1,sample_frac=1,ref_method = 'median'))
})

test_that("PQN throws an error if sample filtering is too strict.",{
    expect_error(out <- pqn_normalisation(df=testData$data, classes=testData$class, 
        qc_label="QC",qc_frac = 0.1,sample_frac=1.1,ref_method = 'median'))
})




