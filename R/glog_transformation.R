#' @importFrom stats optimize
#'

NULL

#' variance stabilising (extended) generalised logarithm transformation
#'
#' https://doi.org/10.1186/1471-2105-8-234
#'
#' @param y values to be tranformed
#' @param y0 offset applied to y (default=0).
#' @param lambda transform parameter
#' @return vector of transformed values

glog <- function(y, y0=0, lambda){
    z <- log((y - y0) + sqrt((y - y0)^2 + lambda))
}

#' Internal function for max. likelihood optimisation of glog params
#' Calculates the alternative Jacobian function described in
#' https://doi.org/10.1186/1471-2105-8-234
#'
#' @param y values.
#' @param y0 offset applied to y (default=0).
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

#' Performs glog transformation on the data set,
#' using QC samples to estimate lambda.
#'
#' https://doi.org/10.1186/1471-2105-8-234
#'
#' @param df Peak intensity matrix
#' @param classes Vector of class labels
#' @param qc_label Class label for QC sample
#' @param store_lambda If value of optimised lambda parameter needs to be
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
    step_threshold <- 1e-16 #stop optimisation when improvement less than this

    VF <- apply(df_qc, 1, var, na.rm=TRUE) # variance of all features
    
    # Use the samllest non-0 variance value as optimisation minimum
    newminVar <- sort(VF)[sort(VF) > 0][1]

    low_lim <- 0 ## set lower  and upper limit of optimisation
    upper_lim <- max(pmax(VF, max(VF) / newminVar))

    # search for optimal value of lambda.
    # NB y0 set to default of 0 as not being implemented here
    lambda <- optimize(f=SSE, interval=c(low_lim, upper_lim), y0=0,
        y=df_qc, tol=step_threshold)

    lambda <- as.numeric(lambda[[1]]) # make sure value of objective is numeric
    lambda_opt <- lambda

    # If SSE optimisation fails use fixed lambda value
    lambda_std <- 5.0278 * 10^(-9)
    
    # if optimisation reached upper limit then trigger use of fixed value
    if (abs(upper_lim - lambda) <= 1e-05) {
        cat("Error!Lambda tending to infinity!Using standard\n")
        error_flag <- TRUE
        # if optimisation reached lower limit then trigger use of fixed value
    } else if (abs(low_lim - lambda) <= 1e-05) {
        cat("Error!Lambda tending to -infinity!Using standard\n")
        error_flag <- TRUE
    } else {
        error_flag <- FALSE
    }

    # if flag triggered then apply scale factor
    if (error_flag) {
        lambda <- lambda_std # set lambda to default
        # set scaling factor to 1 / mean (total signal) over all samples
        scal_fact <- apply(df, 2, sum, na.rm=TRUE)
        scal_fact <- mean(scal_fact)
        scal_fact <- 1 / scal_fact
        df <- df * scal_fact # apply scaling factor
    }

    df <- df - min(df, na.rm=TRUE) # set minimum over all values to 0
    df_glog <- as.data.frame(glog(df, 0, lambda)) # apply glog using optimised

    if (store_lambda){
        return(list(df_glog, lambda, lambda_opt, error_flag))
    } else {
        return(df_glog)
    }
}
