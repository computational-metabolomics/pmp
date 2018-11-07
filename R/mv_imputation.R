#' @importFrom  impute impute.knn
#' @importFrom missForest missForest
#' @importFrom pcaMethods pca
NULL

#' Missing value imputation using different algorithms
#'
#' @param df A peak matrix with features in the rows, samples in the columns
#' @param method Missing value imputation method. Supported methods are "knn", "rf", "bpca", "sv", "mn" and "md".
#' @param k Number of neighbors to be used in the imputation 
#' @param rowmax Fraction of missing values per row.
#' @param colmax Fraction of missing values per column.
#' @param maxp Number of features to run on single core. If set to NULL will use total number of features.
#' @param check_df If set to TRUE will check if input data needs to be transposed, so that features are in rows.
#' @export

mv_imputation = function(df, method, k=10, rowmax=0.5, colmax=0.5, maxp=NULL, check_df=TRUE){
  
  if (check_df ==T){
    
    df <- check_peak_matrix_orientation(peak_data = df)
  }
  
  if (is.null(maxp)){
    
    maxp <- max(dim(df))
  }

  if(any(apply(df, 1, function(vec) all(is.na(vec))) == TRUE)) {
    stop("Error occurred. Rows with 100% missing values detected - please remove these rows using the peak filter tool")
  }

  if(any(apply(df, 2, function(vec) all(is.na(vec))) == TRUE)) {
    stop("Error occurred. Columns with 100% missing values detected - please remove these columns using the sample filter tool")
  }

  if (tolower(method) == "knn"){
  
    obj <- suppressWarnings(impute.knn(as.matrix(df), k=k, rowmax=rowmax, colmax=colmax,maxp = maxp))
    df <- obj$data
    
  } else if (tolower(method) == "rf"){
    
    mf_out = missForest(t(df))
    print(mf_out$OOBerror)
    df = t(mf_out$ximp)
    
  } else if (tolower(method) == "bpca"){
    
    pcaOb = pcaMethods::pca(t(df), method ="bpca", scale="none")
    df = t(pcaOb@completeObs)
    df[df<0] = min(df [df>0]) ##GUARD AGAINST NEGATIVE VALUES
    
  } else if (tolower(method) == "sv"){
    
    df[is.na(df)] = min(df, na.rm=TRUE)/2 ##SMV
    
  } else if (tolower(method) == "mn"){
    
    meanrep = function(mat) apply(mat, 1, mean, na.rm=TRUE) ###MEAN REP
    meanVec = meanrep(df)
    for (i in 1:(nrow(df))){
      df[i, ][is.na(df[i, ])] = meanVec[i]
    }
    
  } else if (tolower(method) == "md"){
    
    medianrep = function(mat) apply(mat, 1, median, na.rm=TRUE)
    medianVec = medianrep(df)
    for (i in 1:(nrow(df))){
      df[i, ][is.na(df[i, ])] = medianVec[i]
    }
    
  } else {
    
    stop("Error occurred. No method selected")
  }
  
  df <- as.data.frame(df)
  return(df)
}
