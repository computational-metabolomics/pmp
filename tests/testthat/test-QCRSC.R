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
  expect_equal (out, testData$QCRSC)
})

test_that ("QC-RSC returns error if QC sample label is wrong", {
  classes[classes == "QC"] <- "Quality"
  expect_error(QCRSC(df = data, order = order, batch = batch, classes = classes,
               spar = 0, minQC = 4))
})
