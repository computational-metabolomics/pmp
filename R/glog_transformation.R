#' @importFrom stats optimise
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_vline
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 theme_bw
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

#' Plot SSE error of lamba optimisation process
#' 
#' @param optimised_lambda optimised lamba value from 'glog_optimise_output'
#' @param data_qc peak intenisty matrix of QC samples
#' @param upper_lim limit of the upper lambda value
#' @return ggplot object containing optimisation plot
#' 

glog_plot_optimised_labmda <- function(optimised_lambda=NA, data_qc, upper_lim){
    sse_df <- data.frame(lambda=seq(0, upper_lim, length.out=100))
    
    k=1
    for (lambda in sse_df[, 1]) {
        sse_df[k, 2] <- SSE(lambda=lambda, y0=0, y=t(data_qc))
        k=k+1
    }
    
    colnames(sse_df) <- c('lambda','SSE')
    
    g=ggplot(data=sse_df, aes(x=lambda,y=SSE))+
        geom_vline(xintercept=optimised_lambda, color="red")+
        geom_line() + theme_bw()
    
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
#' @param store_lambda if value of optimised lambda parameter needs to be
#'returned
#' @examples
#' attach (testData)
#' out <- mv_imputation(df=testData$data, method='knn')
#' out <- glog_transformation (df=out, classes=testData$class,
#'     qc_label='QC')
#'
#' @return data frame, peak inntensity matrix after glog transformation
#' @export glog_transformation

glog_transformation <- function(df, classes, qc_label, store_lambda=FALSE) {
    # check if qc_label is present in the classes vector
    if (length(which(classes %in% qc_label)) == 0) {
        stop("QC sample label is not present. Check your qc_label parameter.")
    }
    df <- check_peak_matrix(peak_data=df, classes=classes)

    df_qc <- df[, classes == qc_label]

    offset <- min(df_qc, na.rm=TRUE) # set offset to the minimum QC samples
    df_qc <- df_qc - offset # set minimum of qc data to 0

    VF <- apply(df_qc, 1, var, na.rm=TRUE) # variance of all features
    
    # Upper limit max var or largest ratio max(var)/min(var)
    upper_lim <- max(pmax(VF, max(VF) / sort(VF)[sort(VF) > 0][1]))

    # search for optimal value of lambda.
    lambda <- glog_omptimise_lambda (upper_lim=upper_lim, df_qc=df_qc)
    lambda <- lambda$minimum
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
        df <- glog_rescale_data(df)
    }

    df <- df - min(df, na.rm=TRUE) # set minimum over all values to 0
    df_glog <- as.data.frame(glog(df, 0, lambda)) # apply glog

    if (store_lambda){
        return(list(data=df_glog, lambda=lambda, lambda_opt=lambda_opt, 
            error_flag=error_flag))
    } else {
        return(df_glog)
    }
}
