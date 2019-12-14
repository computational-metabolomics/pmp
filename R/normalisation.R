#' Normalise peak table to the total sum of peak intensities
#' @param df data frame
#' @param check_df ff set to TRUE will check if input data needs to be
#'transposed, so that features are in rows
#' @return normalised peak matrix
#' @examples 
#' out <- normalise_to_sum (pmp:::testData$data)
#'
#' @export

normalise_to_sum <- function(df, check_df=TRUE) {
    if (check_df == TRUE) {
        df <- check_peak_matrix(peak_data=df)
    }
    return(sweep(df, 2, colSums(df, na.rm=TRUE)/100, FUN="/"))
}


#' Calculate reference mean of  samples
#' 
#' @param df_qc peak matrix of QC samples
#' 
#' @return vector of reference mean values
calculate_ref_mean <- function(df_qc){
    ref_mean <-apply(df_qc, 1, mean, na.rm=TRUE)
    return(ref_mean)
}

#' Normalise peak table using PQN method
#'
#' @param df data frame
#' @param classes vector of class labels
#' @param qc_label label used for QC samples. If set to 'all', all samples will
#'be used to calculate correction factor
#' @param ref_mean Vector of reference mean values to use instead of calculating
#' from QC sample group. If set to NULL, QC sample data will be used.
#' @return list of normalised data set and correction coefficients
#' @examples 
#' attach (pmp:::testData)
#' pqn_normalisation(df=pmp:::testData$data,
#'     classes=pmp:::testData$class, qc_label='QC')
#' 
#' @export

pqn_normalisation <- function(df, classes, qc_label, ref_mean=NULL) {
    
    df <- check_peak_matrix(peak_data=df, classes=classes)
    
    if (is.null(ref_mean)){
        if (qc_label == "all") {
            ref <- df
        } else {
            ref <- df[, classes == qc_label]
        }
        ref_mean <- calculate_ref_mean(df_qc=ref)
    }
    
    coef <- vector()
    
    for (i in seq_len(dim(df)[2])) {
        tempMat <- cbind(ref_mean, df[, i])
        vecelim <- which(apply(tempMat, 1, function(x) any(is.na(x))))
        
        if (length(vecelim) != 0) {
            tempMat <- tempMat[-c(vecelim), , drop=FALSE]
        }
        
        coef[i] <- median(as.numeric(tempMat[, 2]/tempMat[, 1]), na.rm=TRUE)
    }
    out <- list(df=df/coef[col(df)], coef=coef)
    return(out)
}
