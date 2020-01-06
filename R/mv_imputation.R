#' @importFrom  impute impute.knn
#' @importFrom missForest missForest
#' @importFrom pcaMethods pca
NULL

#' Replace missing values for specific peak matrix feature
#' @param x index of matrix feature
#' @param df a peak matrix with features in the rows, samples in the columns
#' @param vals vector object with the same length as number of features in 
#'peak matrix to replca NA values with
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
#' @param df a peak matrix with features in the rows, samples in the columns
#' @param method mode to use for missing value imputation. 'mn' for mean and 
#''md'for median.
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
#' @param df a peak matrix with features in the rows, samples in the columns
#' @param method missing value imputation method. Supported methods are 'knn',
#''rf', 'bpca', 'sv', 'mn' and 'md'
#' @param k number of neighbors to be used in the imputation 
#' @param rowmax fraction of missing values per row
#' @param colmax fraction of missing values per column
#' @param maxp number of features to run on single core. If set to NULL will
#'use total number of features
#' @param check_df if set to TRUE will check if input data needs to be
#'transposed, so that features are in rows
#' @return data frame of missing value imputed peak intensity matrix
#' 
#' @examples 
#' df <- MTBLS79[ , MTBLS79$Batch==1]
#' out <- mv_imputation(df=df, method='knn')
#' 
#' @export

mv_imputation <- function(df, method, k=10, rowmax=0.5, colmax=0.5, 
    maxp=NULL, check_df=TRUE) {
    if (check_df == TRUE) {
        df <- check_input_data(df=df)
    }
    if (is.null(maxp)) {
        maxp <- max(dim(df))
    }
    if (any(apply(assay(df), 1, function(vec) all(is.na(vec))) == TRUE)) {
        stop("Error occurred. Rows with 100% missing values detected - please
    remove these rows using the peak filter tool")
    }
    if (any(apply(assay(df), 2, function(vec) all(is.na(vec))) == TRUE)) {
        stop("Error occurred. Columns with 100% missing values detected - please
    remove these columns using the sample filter tool")
    }
    if (tolower(method) == "knn") {
        obj <- suppressWarnings(impute.knn(assay(df), k=k, 
            rowmax=rowmax, colmax=colmax, maxp=maxp))
        assay(df) <- obj$data
    } else if (tolower(method) == "rf") {
        mf_out <- missForest(t(assay(df)))
        print(mf_out$OOBerror)
        assay(df) <- t(mf_out$ximp)
    } else if (tolower(method) == "bpca") {
        pcaOb <- pcaMethods::pca(t(assay(df)), method="bpca", scale="none")
        bpca_out <- t(pcaOb@completeObs)
        ##GUARD AGAINST NEGATIVE VALUES
        bpca_out[bpca_out < 0] <- min(bpca_out[bpca_out > 0])  
        assay(df) <- bpca_out
    } else if (tolower(method) == "sv") {
        assay(df)[is.na(assay(df))] <- min(assay(df), na.rm=TRUE)/2  ##SMV
    } else if (tolower(method) == "mn") {
        assay(df) <- impute_mode(df=assay(df), method="mn")
    } else if (tolower(method)=="md") {
        assay(df) <- impute_mode(df=assay(df), method="md")
    } else {
        stop("Error occurred. No method selected")
    }
    meta_data <- metadata(df)
    meta_data$processing_history$mv_imputation <- return_function_args()
    metadata(df) <- meta_data
    df <- return_original_data_structure(df)
    return(df)
}
