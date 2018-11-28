context("test-dorsdtable")

test_that("RSD table is returned in proper order", {
  out <- doRSDtable(testData$doRSD)
  expect_equal(out , testData$doRSDtable)
})

test_that("RSD table is returned when NA's are not present in data", {
  data <- testData$data
  data[is.na(data)] <- 0
  # In boht blank samples this feature has NA, so will couse NaN
  data <- data[-5,]
  rsd <- doRSD(Data=data, classes = testData$class)
  out <- doRSDtable(rsd)
  expect_equal(out , testData$doRSDtable_no_NA)
})