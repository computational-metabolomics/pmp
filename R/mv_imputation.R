#' @importFrom  impute impute.knn
#' @importFrom missForest missForest
#' @importFrom pcaMethods pca
#' @importFrom matrixStats rowMedians
NULL

#' Replace missing values for specific peak matrix feature.
#' @param x \code{integer(1)}, index of matrix feature.
#' @param df a \code{numeric()}, peak matrix with features in the rows,
#' samples in the columns.
#' @param vals \code{numeric()}, vector object with the same length as number
#' of features in peak matrix to replca NA values with.
#' 
#' @return feature at specified index 'x' with replaced missing values
#' 
#' @noRd
replace_na <- function(x, df, vals){
    col_vals <- df[x, ]
    na_val <- vals[x]
    col_vals[is.na(col_vals)] <- na_val
    col_vals
}

#' Imput missing values using mode of each feature
#' @param df \code{numeric()}, a peak matrix with features in the rows, samples
#' in the columns
#' @param method \code{character(1)}, mode to use for missing value imputation.
#' 'mn' for mean and md'for median.
#' 
#' @return data frame of missing value imputed peak intensity matrix
#' 
#' @noRd
impute_mode <- function (df, method){
    if (method=="mn"){
        replacement_vals <- rowMeans(df, na.rm=TRUE)
    }
    if (method=="md"){
        replacement_vals <- rowMedians(df, na.rm=TRUE)
    }
    feature_names <- rownames(df)
    df <- do.call(rbind, lapply(seq_len(nrow(df)), replace_na, df=df, 
        vals=replacement_vals))
    rownames(df) <- feature_names
    df
}

#' Missing value imputation using different algorithms
#' 
#' Missing values in metabolomics data sets occur widely and can originate from 
#' a number of sources, including technical and biological reasons. \cr
#' Missing values imputation is applied to replace non-existing values 
#' with an estimated values while maintaining the data structure. A number of 
#' different methods are available as part of this function. \cr
#' 
#' Supported missing value imputation methods are: \cr
#' \cr
#' \code{knn} - K-nearest neighbour. For each feature in each sample, missing 
#' values are replaced by the mean average value (non-weighted) calculated 
#' from its \code{k} closest neighbours in multivariate space (default distance 
#' metric: euclidean distance); \cr
#' \cr
#' \code{rf} - Random Forest. This method is a wrapper of 
#' \link[missForest]{missForest} function. For each feature, missing values are 
#' iteratively imputed until a maximum number of iterations (10), or until the
#' difference between consecutively-imputed matrices becomes positive. 
#' Trees per forest are set to 100, variables included per tree are calculate 
#' using formula \eqn{sqrt(total number of variables)}; \cr
#' \cr
#' \code{bpca} - Bayesian principal component analysis. This method is a 
#' wrapper of \link[pcaMethods]{pca} function. Missing values are replaced by
#' the values obtained from principal component analysis regression with a 
#' Bayesian method. Therefore every imputed missing value does not occur 
#' multiple times, neither across the samples nor across the metabolite 
#' features; \cr
#' \cr
#' \code{sv} - Small value. For each feature, replace missing values with half
#' of the lowest value recorded in the entire data matrix; \cr
#' \cr
#' \code{'mn'} - Mean. For each feature, replace missing values with the mean 
#' average (non-weighted) of all other non-missing values for that variable;\cr
#' \cr
#' \code{'md'} - Median. For each feature, replace missing values with the 
#' median of all other non-missing values for that variable. \cr
#' \cr
#' 
#' @inheritParams filter_peaks_by_blank
#' @param method \code{character(1)}, missing value imputation method.
#' Supported methods are \code{knn}, \code{rf}, \code{bpca}, \code{sv},
#' \code{'mn'} and \code{'md'}.
#' @param k \code{numeric(1)}, for a given sample containing a missing value,
#' the number of nearest neighbours to include to calculate a replacement 
#' value. Used only for method \code{knn}. 
#' @param rowmax \code{numeric(1)}, the maximum percentage of missing data 
#' allowed in any row. For any rows exceeding given limit, missing values are 
#' imputed using the overall mean per sample. Used only for method \code{knn}.
#' @param colmax \code{numeric(1)}, the maximum percent missing data allowed in
#' any column. If any column exceeds given limit, the function will report an
#' error Used only for method \code{knn}.
#' @param maxp \code{integer(1)}, number of features to run on single core.
#' If set to NULL will use total number of features.
#' @param check_df \code{logical(1)}, if set to TRUE will check if input data 
#' needs to be transposed, so that features are in rows.
#' @return Object of class \code{SummarizedExperiment}. If input data are a 
#' matrix-like (e.g. an ordinary matrix, a data frame) object, function returns 
#' the same R data structure as input with all value of data type 
#' \code{numeric()}.
#' 
#' @examples 
#' df <- MTBLS79[ ,MTBLS79$Batch == 1]
#' out <- mv_imputation(df=df, method='knn')
#' 
#' @export

mv_imputation <- function(df, method, k=10, rowmax=0.5, colmax=0.5, 
    maxp=NULL, check_df=TRUE) {
    if (check_df == TRUE) {
        df <- check_input_data(df=df)
    } else {
        # Avoid transposing peak matrix using fake class label
        df <- check_input_data(df=df, classes=rep(1, ncol(df)))
    }
    if (is.null(maxp)) {
        maxp <- max(dim(df))
    }
    if (any(rowSums(is.na(assay(df))) == ncol(df))) {
        stop("Error occurred. Rows with 100% missing values detected - please
    remove these rows using the peak filter tool")
    }
    if (any(colSums(is.na(assay(df))) == nrow(df))) {
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
