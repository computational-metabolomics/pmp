context ("Test QCRSC function")

classes <- MTBLS79$Class
batch <- MTBLS79$Batch
order <- c(1:ncol(MTBLS79))
features_to_select <- c("70.03364", "70.03375", "70.03413", "70.03427", 
  "73.0584", "73.05882", "73.53802", "73.53822", "73.53838", "75.48804",
  "75.48817", "75.48845", "76.8677", "76.86788", "76.86806", "78.02111", 
  "97.00463", "97.005", "98.99486", "98.99512")
data <- assay(MTBLS79[features_to_select, ])

test_that ("QC-RSC returns expected output", {
  out <- QCRSC(df = data, order = order, batch = batch, classes = classes,
                  spar = 0, minQC = 4)
  attributes(out)$processing_history <- NULL
  expect_equal (out, testData$QCRSC, tolerance = 10^-5)
})

test_that ("QC-RSC returns error if QC sample label is wrong.", {
  classes[classes == "QC"] <- "Quality"
  expect_error(QCRSC(df = data, order = order, batch = batch, classes = classes,
               spar = 0, minQC = 4))
})

test_that ("QC-RSC returns the expected reference when batch_ref = 'median_all'.", {
    out <- QCRSC(df = data, order = order, batch = batch, classes = classes,
        spar = 0, minQC = 4,batch_ref = 'median_all')
    A=attributes(out)
    expect_equal(A$processing_history$QCRSC$computed_batch_ref[1],median(data[1,],na.rm=TRUE),tolerance = 0.0005)
})

test_that ("QC-RSC returns the expected reference when batch_ref = 'median_qc'.", {
    out <- QCRSC(df = data, order = order, batch = batch, classes = classes,
        spar = 0, minQC = 4,batch_ref = 'median_qc')
    A=attributes(out)
    
    expect_equal(A$processing_history$QCRSC$computed_batch_ref[1],median(data[1,classes=='QC'],na.rm=TRUE),tolerance = 0.0005)
})

test_that ("QC-RSC returns the expected reference when batch_ref = 'median_sample'.", {
    out <- QCRSC(df = data, order = order, batch = batch, classes = classes,
        spar = 0, minQC = 4,batch_ref = 'median_sample')
    A=attributes(out)
    expect_equal(A$processing_history$QCRSC$computed_batch_ref[1],median(data[1,classes!='QC'],na.rm=TRUE),tolerance = 0.0005)
})

test_that ("QC-RSC returns the expected reference when batch_ref = 100.", {
    out <- QCRSC(df = data, order = order, batch = batch, classes = classes,
        spar = 0, minQC = 4,batch_ref = 100)
    A=attributes(out)
    expect_true(all(A$processing_history$QCRSC$computed_batch_ref==100))
    expect_true(length(A$processing_history$QCRSC$computed_batch_ref)==20)
})

test_that ("QC-RSC returns the expected reference when batch_ref = c(1:20).", {
    out <- QCRSC(df = data, order = order, batch = batch, classes = classes,
        spar = 0, minQC = 4,batch_ref = c(1:20))
    A=attributes(out)
    expect_true(all(A$processing_history$QCRSC$computed_batch_ref==1:20))
    expect_true(length(A$processing_history$QCRSC$computed_batch_ref)==20)
})

test_that ("QC-RSC throws an error if the reference is the wrong length.", {
    expect_error(QCRSC(df = data, order = order, batch = batch, classes = classes,
        spar = 0, minQC = 4,batch_ref = c(100,100)))
})

test_that ("QC-RSC throws an error if batch_ref is not a suitable string.", {
    expect_error(QCRSC(df = data, order = order, batch = batch, classes = classes,
        spar = 0, minQC = 4,batch_ref = 'cake'))
})

test_that ("QC-RSC throws an error if batch_method is not a suitable string.", {
    expect_error(QCRSC(df = data, order = order, batch = batch, classes = classes,
        spar = 0, minQC = 4,batch_ref = 100, batch_method='cake'))
})


test_that ("QC-RSC returns expected values when using batch_method = 'offset'.", {
    out=QCRSC(df = data, order = order, batch = batch, classes = classes,
        spar = 0, minQC = 4,batch_ref = 100, batch_method='offset')
    expect_equal(out[1,1],-35.40889,tolerance=5e-6)
    expect_equal(out[2,2],126.6152,tolerance=5e-6)
})

test_that ("QC-RSC replaces 0 with expected values", {
    
    # 0 in first sample
    data[1,1]=0
    
    # replace 0 with NA using QCRMS
    out=QCRSC(df = data, order = order, batch = batch, classes = classes,
        spar = 0, minQC = 4,batch_ref = 'median_all', batch_method='offset',replace_zero = NA)
    
    A=attributes(out)
    
    # QCRMS should have put the 0 back when finished
    expect_equal(out[1,1],0)
    
    # QCRMS should have replaced 0 with NA, so compare QCRMS with manually replaced
    data[1,1]=NA
    out=QCRSC(df = data, order = order, batch = batch, classes = classes,
        spar = 0, minQC = 4,batch_ref = 'median_all', batch_method='offset')
    
    B = attributes(out)
    expect_equal(
        A$processing_history$QCRSC$computed_batch_ref[1],
        B$processing_history$QCRSC$computed_batch_ref[1],tolerance = 0.0005)
    
    
    # 0 in first sample
    data[1,1]=0
    
    # replace 0 with 100 using QCRMS
    out=QCRSC(df = data, order = order, batch = batch, classes = classes,
        spar = 0, minQC = 4,batch_ref = 'median_all', batch_method='offset',replace_zero = 100)
    
    A=attributes(out)
    
    # QCRMS should have put the 0 back when finished
    expect_equal(out[1,1],0)
    
    # QCRMS should have replaced 0 with NA, so compare QCRMS with manually replaced
    data[1,1]=100
    out=QCRSC(df = data, order = order, batch = batch, classes = classes,
        spar = 0, minQC = 4,batch_ref = 'median_all', batch_method='offset')
    
    B = attributes(out)
    expect_equal(
        A$processing_history$QCRSC$computed_batch_ref[1],
        B$processing_history$QCRSC$computed_batch_ref[1],tolerance = 0.0005)
})





