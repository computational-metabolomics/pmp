context("test-dorsdtable")

test_that("RSD table is returned in proper order", {
  out <- doRSDtable(testData$doRSD)
  expect_equal(out , testData$doRSDtable)
})
