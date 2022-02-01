context ("sbcWrapper function")

classes <- MTBLS79$Class
batch <- MTBLS79$Batch
order <- c(1 : ncol(MTBLS79))

Data <- assay(MTBLS79)
qcData <- Data[,classes == "QC"]
qc_batch <- batch[classes == "QC"]
qc_order <- order[classes == "QC"]

test_that ("sbcWrapper returns expected output", {

  out <- sbcWrapper(id=4, qcData=qcData, order=order, qcOrder=qc_order, qcBatch=qc_batch,
                       log=TRUE, spar=0, batch=batch, minQC=4)

  expect_equal (out, testData$sbcWrapper, tolerance=10^-7)

})

test_that ("sbcWrapper returns expected output, with predefined spar value", {
  
  out <- sbcWrapper(id=4, qcData=qcData, order=order, qcOrder=qc_order, qcBatch=qc_batch,
                    log=TRUE, spar=1, batch=batch, minQC=4)
  
  expect_equal (out, testData$sbcWrapper_spar_1)
  
})

test_that ("sbcWrapper returns NA if not enough data points are measured", {
  
  qcData[4, -c(sample(38, 3))] <- NA
  
  out <- sbcWrapper(id=4, qcData=qcData, order=order, qcOrder=qc_order, qcBatch=qc_batch,
                    log=TRUE, spar=1, batch=batch, minQC=4)
  
  expect_true (all(is.na(out)))
  
})

test_that ("sbcWrapper returns error", {
  dim(qcData) <- NULL
  expect_error(sbcWrapper(id=4, qcData=qcData, order=order, qcOrder=qc_order, qcBatch=qc_batch,
                    log=TRUE, spar=1, batch=batch, minQC=4))
})



