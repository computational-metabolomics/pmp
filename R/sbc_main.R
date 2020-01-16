#' @importFrom matrixStats rowMedians
NULL

#' Quality Control-Robust Spline Correction (QC-RSC)
#' 
#' Implementation of Quality QC-RSC algorithm for signal drift and batch 
#' effect correction within/across a multi-batch direct infusion mass
#' spectrometry (DIMS) and liquid chromatography mass spectrometry (LCMS)
#' datasets.
#' This version supports missing values, but requires at least 4 data point 
#' for quality control (QC) samples measured within each analytical batch.
#' The smoothing parameter (spar) can be optimised using leave-one-out cross
#' validation to avoid overfitting.
#' 
#' @author Andris Jankevics \email{a.jankevics@bham.ac.uk}
#' @references Kirwan et al, Anal. Bioanal. Chem., 405 (15), 2013
#'\url{https://dx.doi.org/10.1007/s00216-013-6856-7}
#' 
#' @inheritParams sbc_plot
#' @param order \code{numeric()}, A numeric vector indicating the order in 
#' which samples were measured.
#' @param spar \code{numeric(1)}, Spline smoothing parameter. Should be in the 
#' range 0 to 1. If set to 0 it will be estimated using leave-one-out 
#' cross-validation.
#' @param log \code{logical(1)}, to perform the signal correction fit on the
#' log scaled data. Default is TRUE.
#' @param minQC \code{integer(1)}, Minimum number of measured quality control
#' (QC) samples required for signal correction within feature per batch. For 
#' features where signal was measured in less QC samples than threshold signal
#' correction won't be applied.
#'
#' @return Object of class \code{SummarizedExperiment}. If input data are a 
#' matrix-like (e.g. an ordinary matrix, a data frame) object, function returns 
#' the same R data structure as input with all value of data type 
#' \code{numeric()}.
#'
#' @examples 
#' classes <- MTBLS79$Class
#' batch <- MTBLS79$Batch
#' order <- c(1:ncol(MTBLS79))
#' 
#' out <- QCRSC(df = MTBLS79[1:10, ], order = order, batch = MTBLS79$Batch,
#' classes = MTBLS79$Class, spar = 0, minQC = 4)
#' 
#' @export

QCRSC <- function(df, order, batch, classes, spar = 0, log = TRUE,
    minQC = 5, qc_label="QC") {
    
    df <- check_input_data(df=df, classes=classes)
    
    if (length(which(classes == qc_label)) <= 0) {
        message("QC samples are not defined! Please see help page for
            parameter qc_label.")
        stop()
    }
    
    message("The number of NA and <= 0 values in peaksData before QC-RSC: ",
        sum(is.na(assay(df)) | assay(df) <= 0))
    
    qcData <- df[, classes == qc_label]
    qc_batch <- batch[classes == qc_label]
    qc_order <- order[classes == qc_label]
    QC_fit <- lapply(seq_len(nrow(df)), sbcWrapper, qcData = assay(qcData), 
        order = order, qcBatch = qc_batch, qcOrder = qc_order, 
        log = log, spar = spar, batch = batch, minQC = minQC)
    QC_fit <- do.call(rbind, QC_fit)
    
    # Median value for each fature, and divide it by predicted value
    mpa <- matrixStats::rowMedians(assay(df), na.rm=TRUE)
    QC_fit <- QC_fit/mpa
    
    # Divide measured value by correction factor
    assay(df) <- assay(df)/QC_fit
    assay(df)[assay(df) <= 0] <- NA
    
    meta_data <- metadata(df)
    meta_data$processing_history$QCRSC <- return_function_args()
    metadata(df) <- meta_data
    df <- return_original_data_structure(df)
    
    return (df)
}
