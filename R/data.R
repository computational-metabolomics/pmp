#' Internal data set for testhat unit tests
#'
#' Subset of DIMS data set to test package functions
#'
#' @source https://www.ebi.ac.uk/metabolights/MTBLS79
#' @noRd
"testData"

#' Direct-infusion mass spectrometry (DIMS) data set
#'
#' Data set of 20 biological (cow vs sheep) serum samples that were analysed
#' repeatedly, in 8 batches across 7 days.
#'
#' Code below includes all commands used to generate MTBLS79 object. \cr
#' 
#' \preformatted{
#' library (openxlsx)
#' library (SummarizedExperiment)
#' 
#' download.file(destfile = "MTBLS79.xlsx", mode="wb",
#' url = "ftp://ftp.ebi.ac.uk/pub/databases/metabolights/studies/public/ 
#'    MTBLS79/Dataset07__SFPM.xlsx") 
#' wb <- openxlsx::loadWorkbook(xlsxFile="MTBLS79.xlsx")
#'
#' MTBLS79 <- list()
#'
#' MTBLS79$assay <- openxlsx::readWorkbook(wb, "data", colNames=T, rowNames=T)
#' 
#' # Last row of the peak matrix represent mean intensities across all samples.
#' MTBLS79$assay <- MTBLS79$assay[-c(nrow(MTBLS79$assay)), ]
#'
#' # Transpose peak matrix, so that features are in rows and samples in columns.
#' MTBLS79$assay <- as.matrix(t(MTBLS79$assay))
#'
#' # Missing values in the input data are stored as 0, replace with NA
#' MTBLS79$assay[MTBLS79$assay == 0] <- NA
#' 
#' rownames(MTBLS79$assay) <- round(as.numeric(rownames(MTBLS79$assay)), 5)
#' 
#' MTBLS79$colData <- openxlsx::readWorkbook(wb, "meta", colNames=T, rowNames=F)
#' MTBLS79$colData <- MTBLS79$colData[-c(nrow(MTBLS79$colData)), 1:4]
#' 
#' MTBLS79 <- SummarizedExperiment(assays=list(MTBLS79$assay), 
#' colData=DataFrame(MTBLS79$colData))
#' }
#'
#' @references Kirwan et al, Scientific Data volume 1,
#'Article number: 140012 (2014) https://www.nature.com/articles/sdata201412
#'
#' @format A \link[SummarizedExperiment]{SummarizedExperiment} object. \cr
#' \code{assay(MTBLS79)} Peak intensities of the DIMS data set. 
#' Contains 172 samples and 2488 features. \cr
#' \code{colData(MTBLS79)} Sample meta data containing 4 columns. \cr
#'     \code{Batch} - sample batch name. \cr
#'     \code{Sample_Rep} - sample replicate code. \cr
#'     \code{Class} - sample class labels. \cr
#'     \code{Class2} - alternative sample class labels grouping together 
#'     replicate samples. \cr
# 
#' @source https://www.ebi.ac.uk/metabolights/MTBLS79
#' 
#' @docType data
#' 
"MTBLS79"
