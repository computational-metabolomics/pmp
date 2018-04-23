#' SMU normalisation
#' @param df Data frame.
#' @return Normalised data frame.

smu_normalisation = function(df){
  return (sweep(df, 1, rowSums(df, na.rm=TRUE) / 100, FUN="/"))
}


#' PQN normalisation
#'
#' @param df Data frame.
#' @param classes Vector of class labels.
#' @param qc_label Label used for QC samples.
#' @return List of normalised data set and correction coefficients
#' @export

pqn_normalisation = function(df, classes, qc_label){

  if (length(classes)!=ncol(df)) {
    cat("Length of sample labels doesn't match data matrix dimensions .Input data frame should be in format: features rows, samples columns.", "\n")
    stop()
  }

  ref = df[,classes == qc_label]

  ref_mean = apply(ref, 1, mean, na.rm=TRUE)
  coef = vector()

  for (i in 1:dim(df)[2]){
    tempMat = cbind(ref_mean, df[,i])
    vecelim = which(apply(tempMat, 1, function(x) any(is.na(x))))

    if (length(vecelim)!=0){
      tempMat = tempMat[-c(vecelim),]
    }

    coef[i] = median(as.numeric(tempMat[, 2] / tempMat[, 1]), na.rm=TRUE)
  }
  out=list(df= df / coef[col(df)], coef=coef)
  return (out)
}
