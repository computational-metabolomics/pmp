context("test-createClassAndColors")

test_that("Colors and ordered factors of sample classes is correct", {

  out <- createClassAndColors(class=testData$class)
  
  expect_equal(testData$createClassAndColors$class, out$class)
  expect_equal(testData$createClassAndColors$manual_colors, out$manual_colors)
  
})

test_that("createClassAndColors doesn't reorder ordered factor again", {
  class <- c(1,1,"QC",1,2,2,2,2,3,3,"QC")
  class <- factor (class, levels = c("3", "2", "1","QC"), ordered=T)
  out <- createClassAndColors(class)
  expect_equal(out, testData$createClassAndColors_ordered)
})

test_that("createClassAndColors uses reainbow colors if number of classes exceeds ", {
  class <- c(1,1,2,3,3,3,4,4,5,5,7,7,8,8,10,10,11,11,12,12)
  out <- createClassAndColors(class)
  expect_equal(out, testData$createClassAndColors_rainbow)
})