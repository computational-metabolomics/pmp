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
#' @param spar_lim A 2 element numeric vector containing the min and max
#'values of spar when searching for an optimum. Default \code{spar_lim = c(-1.5,1.5)}
#' @param batch_ref Numeric reference value(s) that all batches are corrected to. Must
#' be of length 1, or length equal to the number of features.
#' Allowed character values include "median_all" (default) to use the median of all input 
#' samples (including QCs) , "median_qc" to use the median of all qc samples,
#' and "median_sample" to use the median of all samples excluding QCs as the
#' reference.
#' @param batch_method The method to use for correcting batches. Allowed values
#' include "ratio" (default) to apply the correction using 
#' division/multiplication or "offset" to apply the correction using 
#' subtraction/addition.
#' @param replace_zero A value used to replace all zero values after correction.
#' Default is NA. Zeros are replaced before fitting, and then replaced at the end.
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
    minQC = 5, qc_label="QC", spar_lim = c(-1.5,1.5), batch_ref='median_all', 
    batch_method='ratio',replace_zero=NA) {
    
    df <- check_input_data(df=df, classes=classes)
    
    found0=which(assay(df) == 0)
    assay(df)[found0] <- replace_zero
    
    if (length(which(classes == qc_label)) <= 0) {
        message("QC samples are not defined! Please see help page for
            parameter qc_label.")
        stop()
    }
    
    qcData <- df[, classes == qc_label]
    qc_batch <- batch[classes == qc_label]
    qc_order <- order[classes == qc_label]
    QC_fit <- lapply(seq_len(nrow(df)), sbcWrapper, qcData = assay(qcData), 
        order = order, qcBatch = qc_batch, qcOrder = qc_order, 
        log = log, spar = spar, batch = batch, minQC = minQC,
        spar_lim = spar_lim)
    QC_fit <- do.call(rbind, QC_fit)
    
    meta_data <- metadata(df)
    meta_data$processing_history$QCRSC <- return_function_args()
    meta_data$processing_history$QCRSC$QC_fit=QC_fit
    
    # Median value for each feature, and divide it by predicted value
    
    if (batch_ref[1]=='median_all') {
        # median of all samples (including QCs)
        mpa <- matrixStats::rowMedians(assay(df), na.rm=TRUE)
    } else if (batch_ref[1]=='median_qc') {
        # median of QC samples
        mpa <- matrixStats::rowMedians(assay(qcData), na.rm=TRUE)
    } else if (batch_ref[1] =='median_sample') {
        sData=df[,classes!=qc_label]
        mpa <- matrixStats::rowMedians(assay(sData), na.rm=TRUE)
    } else if (is.numeric(batch_ref)) {
        # assume mpa is provided, so use it
        if (length(batch_ref)>1) {
            # check length
            if (length(batch_ref) != nrow(df)){
                stop('if numeric length(batch_ref) must be equal to 1, or equal to the number of features.')
            }
            mpa=batch_ref
        } else {
            mpa=rep(batch_ref,nrow(df))
        }
        
    } else {
        stop('Provided batch_ref must be numeric, "median_all", "median_qc" or "median_sample".')
    }
    
    # convert mpa to matrix
    mpa = matrix(mpa,nrow=length(mpa),ncol=ncol(df),byrow=FALSE)
    
    if (batch_method == 'ratio') {
        # Divide measured value by correction factor
        assay(df) <- assay(df)/QC_fit
        # scale up to desired value
        assay(df) <- assay(df)*mpa
        
    } else if (batch_method == 'offset') {
        # subtract fit, add offset
        assay(df) <- (assay(df)-QC_fit) + mpa
    } else {
        stop('batch_method must be "offset" or "ratio".')
    }
    
    if (length(found0)>0) {
        assay(df)[found0] <- 0
    }
    
    meta_data$processing_history$QCRSC$computed_batch_ref=mpa[,1]
    
    metadata(df) <- meta_data
    df <- return_original_data_structure(df)
    
    return (df)
}
