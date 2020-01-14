context ("test-splineSmoother function")

y_qc <- c(69573.34, 70434.15, 71633.42, 79970.91, 91370.66, 89932.12,
          83033.16)
x_qc <- c(1, 2, 3, 7, 13, 18, 23)
x_sample <- c(4, 5, 6, 8, 9, 10, 11, 12, 14, 15, 16, 17, 19, 20, 21, 22)

test_that ("splineSmoother returns expected output", {
  out <- splineSmoother(x=x_qc, y=y_qc, newX=x_sample,
    log=TRUE, spar=0, a=1)
  expect_equal (out, testData$splineSmoother)

})
