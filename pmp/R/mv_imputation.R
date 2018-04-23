#' @importFrom  impute impute.knn
#' @importFrom missForest missForest
#' @importFrom pcaMethods pca
NULL

#' Missing value imputation using different algorithms
#'
#' @param df Data frame.
#' @param method Missing value imputation method.
#' @param k Number of neighbour values to use.
#' @param rowmax Fraction of missing values per row.
#' @param colmax Fraction of missing values per column.
#' @export

mv_imputation = function(df, method, k=10, rowmax=0.5, colmax=0.5){

  if(any(apply(df, 1, function(vec) all(is.na(vec))) == TRUE)) {
    stop("Error occurred. Rows with 100% missing values detected - please remove these rows using the peak filter tool")
  }

  if(any(apply(df, 2, function(vec) all(is.na(vec))) == TRUE)) {
    stop("Error occurred. Columns with 100% missing values detected - please remove these columns using the peak filter tool")
  }

  if (tolower(method) == "knn"){
    #suppressWarnings( suppressPackageStartupMessages( stopifnot( library("impute", quietly=TRUE, logical.return=TRUE, character.only=TRUE))))
    obj = suppressWarnings(impute.knn(as.matrix(t(df)), k=k, rowmax=rowmax, colmax=colmax))
    df = as.data.frame(t(obj$data))
  } else if (tolower(method) == "rf"){
    #suppressWarnings( suppressPackageStartupMessages( stopifnot( library("missForest", quietly=TRUE, logical.return=TRUE, character.only=TRUE))))
    mf_out = missForest(df)
    print(mf_out$OOBerror)
    df = mf_out$ximp
  } else if (tolower(method) == "bpca"){
    #suppressWarnings( suppressPackageStartupMessages( stopifnot( library("pcaMethods", quietly=TRUE, logical.return=TRUE, character.only=TRUE))))
    pcaOb = pcaMethods::pca(df, method ="bpca", scale="none")
    df = pcaOb@completeObs
    df[df<0] = min(df [df>0]) ##GUARD AGAINST NEGATIVE VALUES
  } else if (tolower(method) == "sv"){
    df[is.na(df)] = min(df, na.rm=TRUE)/2 ##SMV
  } else if (tolower(method) == "mn"){
    meanrep = function(mat) apply(mat, 2, mean, na.rm=TRUE) ###MEAN REP
    meanVec = meanrep(df)
    for (i in 1:(ncol(df))){
      df[, i][is.na(df[, i])] = meanVec[i]
    }
  } else if (tolower(method) == "md"){
    medianrep = function(mat) apply(mat, 2, median, na.rm=TRUE)
    medianVec = medianrep(df)
    for (i in 1:(ncol(df))){
      df[, i][is.na(df[, i])] = medianVec[i]
    }
  } else {
    stop("Error occurred. No method selected")
  }
  return(df)
}
