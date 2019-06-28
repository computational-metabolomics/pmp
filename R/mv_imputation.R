#' @importFrom  impute impute.knn
#' @importFrom missForest missForest
#' @importFrom pcaMethods pca
NULL

#' Replace missing values for specific peak matrix feature
#' @param x Index of matrix feature
#' @param df A peak matrix with features in the rows, samples in the columns
#' @param vals Vector object with the same length as number of features in 
#' peak matrix to replca NA values with
#' 
#' @return feature at specified index 'x' with replaced missing values
#' 
replace_na <- function(x, df, vals){
    col_vals <- df[x, ]
    na_val <- vals[x]
    col_vals[is.na(col_vals)] <- na_val
    col_vals
}

#' Imput missing values using mode of each feature
#' @param df A peak matrix with features in the rows, samples in the columns
#' @param method Mode to use for missing value imputation. 'mn' for mean and 
#' 'md'for median.
#' 
#' @return data frame of missing value imputed peak intensity matrix
#' 
impute_mode <- function (df, method){
    if (method=="mn"){
        replacement_vals <- apply(df, 1, mean, na.rm=TRUE)
    }
    
    if (method=="md"){
        replacement_vals <- apply(df, 1, median, na.rm=TRUE)
    }
    feature_names <- rownames(df)
    df <- do.call(rbind, lapply(seq_len(nrow(df)), replace_na, df=df, 
        vals=replacement_vals))
    rownames(df) <- feature_names
    df
}

#' Missing value imputation using different algorithms
#'
#' @param df A peak matrix with features in the rows, samples in the columns
#' @param method Missing value imputation method. Supported methods are 'knn',
#'  'rf', 'bpca', 'sv', 'mn' and 'md'.
#' @param k Number of neighbors to be used in the imputation 
#' @param rowmax Fraction of missing values per row.
#' @param colmax Fraction of missing values per column.
#' @param maxp Number of features to run on single core. If set to NULL will
#'  use total number of features.
#' @param check_df If set to TRUE will check if input data needs to be
#'  transposed, so that features are in rows.
#' @return data frame of missing value imputed peak intensity matrix
#' 
#' @examples 
#' attach (testData)
#' out <- mv_imputation(df=t(testData$data), method='knn')
#' 
#' @export

mv_imputation <- function(df, method, k=10, rowmax=0.5, colmax=0.5, 
    maxp=NULL, check_df=TRUE) {
    
    if (check_df == TRUE) {
        df <- check_peak_matrix(peak_data=df)
    }
    
    if (is.null(maxp)) {
        maxp <- max(dim(df))
    }
    
    if (any(apply(df, 1, function(vec) all(is.na(vec))) == TRUE)) {
        stop("Error occurred. Rows with 100% missing values detected - please
    remove these rows using the peak filter tool")
    }
    
    if (any(apply(df, 2, function(vec) all(is.na(vec))) == TRUE)) {
        stop("Error occurred. Columns with 100% missing values detected - please
    remove these columns using the sample filter tool")
    }
    
    if (tolower(method) == "knn") {
        obj <- suppressWarnings(impute.knn(as.matrix(df), k=k, 
            rowmax=rowmax, colmax=colmax, maxp=maxp))
        df <- obj$data
    } else if (tolower(method) == "rf") {
        mf_out <- missForest(t(df))
        print(mf_out$OOBerror)
        df <- t(mf_out$ximp)
    } else if (tolower(method) == "bpca") {
        pcaOb <- pcaMethods::pca(t(df), method="bpca", scale="none")
        df <- t(pcaOb@completeObs)
        df[df < 0] <- min(df[df > 0])  ##GUARD AGAINST NEGATIVE VALUES
    } else if (tolower(method) == "sv") {
        df[is.na(df)] <- min(df, na.rm=TRUE)/2  ##SMV
    } else if (tolower(method) == "mn") {
        df <- impute_mode(df=df, method="mn")
    } else if (tolower(method)=="md") {
        df <- impute_mode(df=df, method="md")
    } else {
        stop("Error occurred. No method selected")
    }
    
    return(as.data.frame(df))
}
