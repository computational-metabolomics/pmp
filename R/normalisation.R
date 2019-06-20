#' Normalise peak table to the total sum of peak intensities
#' @param df Data frame
#' @param check_df If set to TRUE will check if input data needs to be transposed, so that features are in rows.
#' @return Normalised peak matrix.
#' @examples 
#' attach (testData)
#' out <- normalise_to_sum (testData$data)
#'
#' @export

normalise_to_sum <- function(df, check_df=TRUE){
  if (check_df == T){
    df <- check_peak_matrix_orientation(peak_data=df)
  }
  return (sweep(df, 2, colSums(df, na.rm=TRUE) / 100, FUN="/"))
}


#' Normalise peak table using PQN method
#'
#' @param df Data frame.
#' @param classes Vector of class labels. 
#' @param qc_label Label used for QC samples. If set to "all", all samples will be used to calculate correction factor
#' @return List of normalised data set and correction coefficients
#' @examples 
#' attach (testData)
#' pqn_normalisation(df=testData$data, classes=testData$class, qc_label = "QC")
#' 
#' @export

pqn_normalisation <- function(df, classes, qc_label){

  df <- check_peak_matrix_orientation(peak_data=df, classes=classes)

  if (qc_label=="all"){
    ref=df
  } else {
    ref=df[ ,classes == qc_label]
  }
    
  ref_mean <- apply(ref, 1, mean, na.rm=TRUE)
  coef <- vector()

  for (i in seq_len(dim(df)[2])){
    tempMat <- cbind(ref_mean, df[,i])
    vecelim <- which(apply(tempMat, 1, function(x) any(is.na(x))))

    if (length(vecelim)!=0){
      tempMat <- tempMat[-c(vecelim),]
    }

    coef[i] <- median(as.numeric(tempMat[ , 2] / tempMat[ , 1]), na.rm=TRUE)
  }
  out <- list(df=df / coef[col(df)], coef=coef)
  return (out)
}