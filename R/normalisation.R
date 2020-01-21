#' @importFrom matrixStats rowAnyMissings

NULL

#' Normalisation by total sum of the features per sample
#' 
#' For each sample, every feature intensity value is divided by the total sum of
#' all feature intensity values measured in that sample (\code{NA} values
#' ignored by default), before multiplication by 100; the unit is \%.
#' 
#' @inheritParams mv_imputation
#' 
#' @return Object of class \code{SummarizedExperiment}. If input data are a 
#' matrix-like (e.g. an ordinary matrix, a data frame) object, function returns 
#' the same R data structure as input with all value of data type 
#' \code{numeric()}.
#' 
#' @examples 
#' df <- MTBLS79[ ,MTBLS79$Batch == 1]
#' out <- normalise_to_sum (df=df)
#'
#' @export

normalise_to_sum <- function(df, check_df=TRUE) {
    if (check_df == TRUE) {
        df <- check_input_data(df)
    } else {
        # normalise_to_sum doesn't need class labels
        # Create generic class label vector to avoid DF to be transposed
        df <- check_input_data(df=df, classes=rep("S", ncol(df)))
    }
    assay(df) <- (sweep(assay(df), 2, colSums(assay(df), na.rm=TRUE)/100, 
        FUN="/"))
    meta_data <- metadata(df)
    meta_data$processing_history$normalise_to_sum <- 
        return_function_args()
        #list (check_df=check_df)
    metadata(df) <- meta_data
    df <- return_original_data_structure(df)
    df
}

#' Calculate reference mean of samples
#' 
#' @param df_qc \code{numeric()}, peak matrix of QC samples.
#' 
#' @return vector of reference mean values
#' @noRd
calculate_ref_mean <- function(df_qc){
    ref_mean <- rowMeans(df_qc, na.rm=TRUE)
    return(ref_mean)
}

#' Probabilistic quotient normalisation (PQN)
#' 
#' For every feature the mean response is calculated across all QC samples. A 
#' reference vector is then generated. The median between the reference vector
#' and every sample is computed obtaining a vector of coefficients related to
#' each sample. Each sample is then divided by the median value of the vector
#' of coefficients; this median value is different for each sample. This 
#' method was adapted by Dieterle et al. (2006) (see references). Its purpose 
#' is to take into account the concentration changes of some metabolite 
#' features that affect limited regions of the data.
#'
#' @references Dieterle F. et al., Anal. Chem., 78(13), 2006. 
#' http://dx.doi.org/10.1021/ac051632c
#'
#' @inheritParams filter_peaks_by_blank 
#' @param ref_mean \code{numeric()} or \code{NULL}, Vector of reference mean
#' values to use instead of calculating from QC sample group. If set to 
#' \code{NULL}, QC sample data will be used.
#' @return Object of class \code{SummarizedExperiment}. If input data are a 
#' matrix-like (e.g. an ordinary matrix, a data frame) object, function returns 
#' the same R data structure as input with all value of data type 
#' \code{numeric()}.
#' 
#' @examples 
#' df <- MTBLS79[ , MTBLS79$Batch==1]
#' pqn_normalisation(df=df,
#'     classes=df$Class, qc_label='QC')
#' 
#' @export

pqn_normalisation <- function(df, classes, qc_label, ref_mean=NULL) {
    df <- check_input_data(df=df, classes=classes)
    if (is.null(ref_mean)){
        if (qc_label == "all") {
            ref <- df
        } else {
            ref <- df[, classes == qc_label]
        }
        ref_mean <- calculate_ref_mean(df_qc=assay(ref))
    }
    coef <- vector()
    for (i in seq_len(dim(df)[2])) {
        tempMat <- cbind(ref_mean, assay(df)[, i])
        vecelim <- which(rowAnyMissings(tempMat))
        if (length(vecelim) != 0) {
            tempMat <- tempMat[-c(vecelim), , drop=FALSE]
        }
        coef[i] <- median(as.numeric(tempMat[, 2]/tempMat[, 1]), na.rm=TRUE)
    }
    assay(df) <- assay(df)/coef[col(assay(df))]
    col_data <- DataFrame(pqn_coef=coef)
    colData(df) <- cbind(colData(df), col_data)
    meta_data <- metadata(df)
    meta_data$processing_history$pqn_normalisation <- return_function_args()
    metadata(df) <- meta_data
    df <- return_original_data_structure(df)
    if (!is(df, "SummarizedExperiment")){
        attributes(df)$flags <- as.matrix(col_data)
    }
    return(df)
}
