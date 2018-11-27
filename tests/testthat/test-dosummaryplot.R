context("test-dosummaryplot")

test_that("dosummaryplot creates ggplot output", {
  out <- doSummaryPlot (Data=testData$data, classes=testData$class, plotTitle="PCA", 
                 blank=NULL, PQN=F, mv_impute=T, glogScaling=T, scale=T, qc_label="QC", 
                 ignorelabel="Removed", labels="all", qc_shape=17, base_size = 10, 
                 pccomp=c(1,2), plot=F)
  expect_equal(out, testData$doSummaryPlot)
})

test_that("dosummaryplot creates ggplot output for several batches", {
  dataList <- list(testData$data[1:15,], testData$data[16:30,])
  classList <- list(testData$class, testData$class)
  titleList <- list("Batch 1", "Batch 2")
  
  out <- doSummaryPlot (Data=dataList, classes=classList, plotTitle=titleList, 
                        blank=NULL, PQN=F, mv_impute=T, glogScaling=T, scale=T, qc_label=NULL, 
                        ignorelabel="Removed", labels="all", qc_shape=17, base_size = 10, 
                        pccomp=c(1,2), plot=F)
  expect_equal(out, testData$doSummaryPlot_multiple_batches)
})

test_that("dosummaryplot stops if input isn't list or data.frame/matrix", {
  expect_error(doSummaryPlot (Data=testData$data[,1], classes=testData$class, plotTitle="PCA", 
                        blank=NULL, PQN=F, mv_impute=T, glogScaling=T, scale=T, qc_label="QC", 
                        ignorelabel="Removed", labels="all", qc_shape=17, base_size = 10, 
                        pccomp=c(1,2), plot=F))
})