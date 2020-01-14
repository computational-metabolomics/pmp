#' @importFrom stats median
NULL

##' @title Quality Control-Robust Spline Correction (QC-RSC)
##' @description Implementation of signal correction algorithm as described
##'in Kirwan et al, Anal. Bioanal. Chem., 405 (15), 2013
##' @author Andris Jankevics \email{a.jankevics@bham.ac.uk}
##' @references Kirwan et al, Anal. Bioanal. Chem., 405 (15), 2013
##'\url{https://dx.doi.org/10.1007/s00216-013-6856-7}
##' @details
##' The smoothing parameter (spar) can be optimised using leave-one-out cross
##'validation to avoid overfitting.
##' 
##' @param df A data frame of values to be corrected (samples in columns and
##'features in rows).
##' @param order A numeric vector indicating the order in which samples
##'were measured.
##' @param batch A vector indicating the batch each sample was measured in.
##'If only one batch then all values should be set to 1.
##' @param classes A factor or character vector of sample classes. All QC
##'samples should be labelled  'QC'.
##' @param spar Spline smoothing parameter. Should be in the range 0 to 1.
##'If set to 0 it will be estimated using leave-one-out cross-validation.
##' @param log TRUE or FALSE to perform the signal correction fit on the log
##'scaled data. Default is TRUE.
##' @param minQC Minimum number of QC samples required for signal correction.
##' @param qc_label Class label for QC sample.
##'
##' @return A data frame of corrected values.
##'
##' @examples 
##' classes <- MTBLS79$Class
##' batch <- MTBLS79$Batch
##' order <- c(1:ncol(MTBLS79))
##' data <- SummarizedExperiment::assay(MTBLS79[1:10, ])
##' 
##' out <- QCRSC(df = data, order = order, batch = batch, classes = classes,
##'spar = 0, minQC = 4)
##' 
##' @export

QCRSC <- function(df, order, batch, classes, spar = 0, log = TRUE,
    minQC = 5, qc_label="QC") {
    
    if (length(which(classes == qc_label)) <= 0) {
        message("QC samples are not defined! Please see help page for
            parameter qc_label.")
        stop()
    }
    
    message("The number of NA and <= 0 values in peaksData before QC-RSC: ",
        sum(is.na(df) | df <= 0))
    
    qcData <- df[, classes == qc_label]
    qc_batch <- batch[classes == qc_label]
    qc_order <- order[classes == qc_label]
    QC_fit <- lapply(seq_len(nrow(df)), sbcWrapper, qcData = qcData, 
        order = order, qcBatch = qc_batch, qcOrder = qc_order, 
        log = log, spar = spar, batch = batch, minQC = minQC)
    QC_fit <- do.call(rbind, QC_fit)
    
    ## Median value for each fature, and divide it by predicted value
    mpa <- apply(df, 1, median, na.rm = TRUE)
    QC_fit <- QC_fit/mpa
    
    ## Divide measured value by correction factor
    res <- df/QC_fit
    res[res <= 0] <- NA
    res
}
