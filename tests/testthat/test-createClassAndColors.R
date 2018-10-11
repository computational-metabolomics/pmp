context("test-createClassAndColors")

test_that("Colors and ordered factors of sample classes is correct", {

  out <- createClassAndColors(class=testData$class)
  
  expect_equal(testData$createClassAndColors$class, out$class)
  expect_equal(testData$createClassAndColors$manual_colors, out$manual_colors)
  
})
