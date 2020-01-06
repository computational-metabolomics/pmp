#' @importFrom stats optimise
#' @import ggplot2
#' @importFrom methods as
#'
NULL

#' Variance stabilising (extended) generalised logarithm transformation
#'
#' https://doi.org/10.1186/1471-2105-8-234
#'
#' @param y values to be tranformed
#' @param y0 offset applied to y (default=0).
#' @param lambda transform parameter
#' @return vector of transformed values

glog <- function(y, y0=0, lambda){
    z <- log((y - y0) + sqrt((y - y0)^2 + lambda))
    return(z)
}

#' Internal function for max. likelihood optimisation of glog params.
#' Calculates the alternative Jacobian function described in
#' https://doi.org/10.1186/1471-2105-8-234
#'
#' @param y values
#' @param y0 offset applied to y (default=0)
#' @param lambda lambda
#' @return numeric, optimised glog parameter

jglog <- function(y, y0=0, lambda){
    z <- glog(y, y0, lambda)
    D <- log(sqrt((y - y0)^2 + lambda))

    # average over all features (bins)
    gmn <- exp(apply(D, 2, mean, na.rm=TRUE))
    zj <- z * gmn  # ML estimate
    return(zj)
}

#' Internal function to estimate SSE for optimising glog params
#'
#' @param y values.
#' @param y0 offset applied to y (default=0)
#' @param lambda transform parameter
#' @return numeric, sum of squared difference

SSE <- function(lambda, y0=0, y) {
    # calculate ML estimate
    z <- jglog(y, y0, lambda)

    # average over all features
    mean_spec <- apply(z, 1, mean, na.rm=TRUE)

    # calculate sum of squared difference between true and estimate
    s <- sum((z - mean_spec)^2, na.rm=TRUE)
    return(s)
}

#' If glog optimisation fails, this function will scale values in the 
#' peak matrix to the 1/mean(total signal) over all samples.
#' 
#' @param df peak intensity matrix
#' @return scaled peak matrix

glog_rescale_data <- function(df){
    # 
    scal_fact <- apply(df, 2, sum, na.rm=TRUE)
    scal_fact <- mean(scal_fact)
    scal_fact <- 1 / scal_fact
    df <- df * scal_fact # apply scaling factor
    return(df)
}

#' Search for optimal value of lambda by minimasing SSE.
#' 
#' @param upper_lim upper limit to use for optimisation
#' @param df_qc peak matrix of QC samples
#'  
#' @return optimised glog lambda value

glog_omptimise_lambda <- function(upper_lim, df_qc){
    lambda <- optimise(f=SSE, interval=c(0, upper_lim), y0=0,
        y=df_qc, tol=1e-16)
    lambda
}

#' Plot SSE error of lambda optimisation process
#' @param df Peak intensity matrix
#' @param optimised_lambda Numeric value of optimised lambda from 
#'     glog_transformation output
#' @param classes vector of class labels
#' @param qc_label class label for QC sample
#' @param plot_grid number of data points to use for SSE optimisation
#' @return ggplot object containing optimisation plot
#' @examples 
#' 
#' classes <- pmp:::testData$class
#' 
#' data <- mv_imputation(df=pmp:::testData$data, method='knn')
#' out <- glog_transformation (df=data, classes=classes,
#'     qc_label='QC')
#' 
#' optimised_lambda <- attributes(out)
#' optimised_lambda <- 
#'     optimised_lambda$processing_history$glog_transformation$lambda_opt
#' 
#' glog_plot_optimised_lambda(df=data, classes=classes,
#'     qc_label="QC", optimised_lambda=optimised_lambda)
#' @export

glog_plot_optimised_lambda <- function(df, optimised_lambda, classes, qc_label,
    plot_grid=100){

    df_qc <- df[, classes == qc_label]
    
    lambda_lim <- c(optimised_lambda, optimised_lambda) +
        c(-optimised_lambda*0.8, optimised_lambda*0.8) 
    
    sse_df <- data.frame(lambda=seq(lambda_lim[1], lambda_lim[2], 
        length.out=plot_grid))
    
    sse_df$SSE <- vapply(X=sse_df$lambda, FUN=SSE, y0=0, y=df_qc, 
        FUN.VALUE=numeric(1))
    
    g <- ggplot(data=sse_df, aes_string(x='lambda',y='SSE')) +
        geom_vline(xintercept=optimised_lambda, color="red") +
        geom_line(size=1.1) + theme_bw() +
        labs (title="Optimisation outup of glog lambda parameter",
            caption=paste("lambda=",optimised_lambda, sep=""))

    return (g)
}

#' Performs glog transformation on the data set,
#' using QC samples to estimate lambda.
#'
#' https://doi.org/10.1186/1471-2105-8-234
#'
#' @param df Peak intensity matrix
#' @param classes vector of class labels
#' @param qc_label class label for QC sample
#' @param lambda If not NULL will use provided value for glog lambda
#' @examples
#' out <- mv_imputation(df=pmp:::testData$data, method='knn')
#' out <- glog_transformation (df=out, classes=pmp:::testData$class,
#'     qc_label='QC')
#'
#' @return An \link[SummarizedExperiment]{SummarizedExperiment} object
#' or numeric object of original input data structure.
#' @export glog_transformation

glog_transformation <- function(df, classes, qc_label, lambda=NULL) {
    # check if qc_label is present in the classes vector
    if (length(which(classes %in% qc_label)) == 0) {
        stop("QC sample label is not present. Check your qc_label parameter.")
    }
    
    df <- check_input_data(df=df, classes=classes)
    df_qc <- df[, classes == qc_label]
    
    offset <- min(assay(df_qc), na.rm=TRUE)#set offset to the minimum QC samples
    assay(df_qc) <- assay(df_qc) - offset # set minimum of qc data to 0
    VF <- apply(assay(df_qc), 1, var, na.rm=TRUE) # variance of all features
    
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

    assay(df) <- assay(df) - min(assay(df), na.rm=TRUE)
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
