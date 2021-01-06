#' @importFrom stats optimise
#' @import ggplot2
#' @importFrom methods as
#' @importFrom matrixStats rowVars
#'
NULL

#' Variance stabilising (extended) generalised logarithm transformation
#'
#' https://doi.org/10.1186/1471-2105-8-234
#'
#' @param y \code{numeric()}, values to be tranformed.
#' @param y0 \code{numeric(1)}, offset applied to y (default=0).
#' @param lambda \code{numeric(1)}, transform parameter.
#' @return vector \code{numeric()}, of transformed values.
#' @noRd

glog <- function(y, y0=0, lambda){
    z <- log((y - y0) + sqrt((y - y0)^2 + lambda))
    return(z)
}

#' Internal function for max. likelihood optimisation of glog params.
#' Calculates the alternative Jacobian function described in
#' https://doi.org/10.1186/1471-2105-8-234
#'
#' @param y \code{numeric()}, values.
#' @param y0 \code{numeric1}, offset applied to y (default=0).
#' @param lambda \code{numeric1}, lambda.
#' @return \code{numeric(1)}, optimised glog parameter.
#' @noRd

jglog <- function(y, y0=0, lambda){
    z <- glog(y, y0, lambda)
    D <- log(sqrt((y - y0)^2 + lambda))
    # average over all features (bins)
    gmn <- exp(colMeans(D, na.rm=TRUE))
    zj <- z * gmn  # ML estimate
    return(zj)
}

#' Internal function to estimate SSE for optimising glog params
#'
#' @param y \code{numeric()}, values.
#' @param y0 \code{numeric(1)}, offset applied to y (default=0).
#' @param lambda \code{numeric(1)}, transform parameter.
#' @return \code{numeric()}, sum of squared difference.
#' @noRd

SSE <- function(lambda, y0=0, y) {
    # calculate ML estimate
    z <- jglog(y, y0, lambda)
    # average over all features
    mean_spec <- rowMeans(z, 1, na.rm=TRUE)
    # calculate sum of squared difference between true and estimate
    s <- sum((z - mean_spec)^2, na.rm=TRUE)
    return(s)
}

#' If glog optimisation fails, this function will scale values in the 
#' peak matrix to the 1/mean(total signal) over all samples.
#' 
#' @param df \code{numeric()}, peak intensity matrix.
#' @return \code{numeric()}, scaled peak matrix.
#' @noRd

glog_rescale_data <- function(df){
    scal_fact <- colSums(df, na.rm=TRUE)
    scal_fact <- mean(scal_fact)
    scal_fact <- 1 / scal_fact
    df <- df * scal_fact # apply scaling factor
    return(df)
}

#' Search for optimal value of lambda by minimasing SSE.
#' 
#' @param upper_lim \code{numeric(1)}, upper limit to use for optimisation.
#' @param df_qc \code{numeric()}, peak matrix of QC samples.
#'  
#' @return \code{numeric(1)}, optimised glog lambda value.
#' @noRd

glog_omptimise_lambda <- function(upper_lim, df_qc){
    lambda <- optimise(f=SSE, interval=c(0, upper_lim), y0=0,
        y=df_qc, tol=1e-16)
    lambda
}

#' Plot SSE error of lambda optimisation process
#' @inheritParams filter_peaks_by_blank
#' @param optimised_lambda \code{numeric(1)}, value of optimised lambda from 
#'  glog_transformation output.
#' @param plot_grid \code{integer(1)}, number of data points to use for SSE 
#' optimisation.
#' @return Class \code{ggplot} object containing optimisation plot.
#' 
#' @examples 
#' data <- MTBLS79[, MTBLS79$Batch == 1]
#' classes <- data$Class
#' 
#' data <- mv_imputation(df=data, method='knn')
#' out <- glog_transformation (df=data, classes=classes,
#'     qc_label='QC')
#' 
#' optimised_lambda <- S4Vectors::metadata(out)
#' optimised_lambda <- 
#'     optimised_lambda$processing_history$glog_transformation$lambda_opt
#' 
#' glog_plot_optimised_lambda(df=data, classes=classes,
#'     qc_label="QC", optimised_lambda=optimised_lambda)
#' @export

glog_plot_optimised_lambda <- function(df, optimised_lambda, classes, qc_label,
    plot_grid=100){
    
    df <- check_input_data(df, classes=classes)
    df_qc <- df[, classes == qc_label]
    
    offset <- min(assay(df_qc), na.rm=TRUE)#set offset to the minimum QC samples
    assay(df_qc) <- assay(df_qc) - offset # set minimum of qc data to 0
    
    lambda_lim <- c(optimised_lambda, optimised_lambda) +
        c(-optimised_lambda*0.8, optimised_lambda*0.8) 
    sse_df <- data.frame(lambda=seq(lambda_lim[1], lambda_lim[2], 
        length.out=plot_grid))
    sse_df$SSE <- vapply(X=sse_df$lambda, FUN=SSE, y0=0, y=assay(df_qc), 
        FUN.VALUE=numeric(1))
    g <- ggplot(data=sse_df, aes_string(x='lambda',y='SSE')) +
        geom_vline(xintercept=optimised_lambda, color="red") +
        geom_line(size=1.1) + theme_bw() +
        labs (title="glog parameter optimisation",
            caption=paste("optimum = ",round(optimised_lambda,2), sep=""))
    return (g)
}

#' Variance stabilising generalised logarithm (glog) transformation
#' 
#' Performs glog transformation on the data set. QC samples can be used to 
#' estimate technical variation in the data set and calculate transformation
#' parameter \eqn{\lambda} (lambda). QC samples usually comprise a pool of 
#' aliquots taken from all other samples in the study  and analysed repeatedly 
#' throughout an analytical batch.
#' 
#' Many univariate and multivariate statistical tests require homogeneity and
#' n ormality of dataset variance. Real-world metabolomics datasets often fail 
#' to meet these criteria due to asymmetric (i.e. non-'normal') and/or 
#' heteroscedatic (i.e. non-homogenous) variance structure. To address this 
#' issue, \code{glog} data transformations may be applied. \cr
#' For each cell within the data matrix, transform the raw value (x) according 
#' to:  \eqn{log10(x + sqrt(x^2 + \lambda))}. The parameter \eqn{\lambda} is 
#' typically calculated using quality control (QC) samples analysed throughout
#' an analysis batch.
#' 
#' @references Parsons HM et. al., BMC Bionf., 8(234), 2007. 
#' https://doi.org/10.1186/1471-2105-8-234
#'
#' @inheritParams filter_peaks_by_blank
#' @param lambda \code{NULL} or \code{numeric(1)}, if not \code{NULL} will use
#' provided value for glog lambda.
#' 
#' @return Object of class \code{SummarizedExperiment}. If input data are a 
#' matrix-like (e.g. an ordinary matrix, a data frame) object, function returns 
#' the same R data structure as input with all value of data type 
#' \code{numeric()}.
#' 
#' @examples
#' df <- MTBLS79[, MTBLS79$Batch == 1]
#' out <- mv_imputation(df=df, method="knn")
#' out <- glog_transformation (df=out, classes=df$Class,
#'     qc_label="QC")
#' 
#' @export glog_transformation

glog_transformation <- function(df, classes, qc_label, lambda=NULL) {
    # check if qc_label is present in the classes vector
    if (length(which(classes %in% qc_label)) == 0) {
        stop("QC sample label is not present. Check your qc_label parameter.")
    }
    df <- check_input_data(df=df, classes=classes)
    df_qc <- df[, classes == qc_label]
    
    # offset for QC samples
    offset <- min(assay(df_qc), na.rm=TRUE)#set offset to the minimum QC samples
    
    assay(df_qc) <- assay(df_qc) - offset # set minimum of qc data to 0
    assay(df) <- assay(df) - offset # apply qc offset to samples
    
    VF <- rowVars(assay(df_qc), na.rm=TRUE) # variance of all features
    # Upper limit max var or largest ratio max(var)/min(var)
    upper_lim <- max(pmax(VF, max(VF) / sort(VF)[sort(VF) > 0][1]))
    if (is.null(lambda)){
        lambda <- glog_omptimise_lambda (upper_lim=upper_lim, 
            df_qc=assay(df_qc))
        lambda <- lambda$minimum
    }
    lambda_opt <- lambda
    error_flag <- FALSE
    # if optimisation reached upper/lower limit then trigger use of fixed value
    if (abs(upper_lim - lambda) <= 1e-05 | abs(0 - lambda) <= 1e-05) {
        cat("Error!Lambda tending to infinity! Using standard\n")
        error_flag <- TRUE
    } 
    # if flag triggered then apply scale factor
    if (error_flag) {
        lambda <- 5.0278 * 10^(-9)
        assay(df) <- glog_rescale_data(assay(df))
    }

    # set minimum over all values to 0
    assay(df) <- glog(assay(df), 0, lambda) # apply glog
    meta_data <- metadata(df)
    meta_data$processing_history$glog_transformation <- 
        c(return_function_args(), list(lambda_opt=lambda_opt, 
            error_flag=error_flag))
    metadata(df) <- meta_data
    df <- return_original_data_structure(df)
    return(df)
}
