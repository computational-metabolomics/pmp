#' variance stabilising generalised logarithm transformation
#'
#' https://doi.org/10.1186/1471-2105-8-234
#'
#' @param y values to be tranformed
#' @param alpha alpha.
#' @param lambda transform parameter
#' @export

glog <- function(y, alpha, lambda){
  z <- log((y - alpha) + sqrt((y - alpha)^2 + lambda))
}

#' Internal glog function
#'
#' @param y values.
#' @param y0 .
#' @param lambda lambda.
#' @export
jglog <- function(y, y0, lambda){
  z <- glog(y, y0, lambda)
  D <- log(sqrt((y - y0)^2 + lambda))

  gmn <- exp(apply(D, 2, mean, na.rm=TRUE))
  zj <- z*gmn
  return(zj)
}

#' Internal glog function to estimate SSE
#'
#' @param y values.
#' @param alpha alpha.
#' @param lambda transform parameter
#' @export

SSE <- function(lambda, alpha, y){
  N <- dim(y)[2]
  len <- dim(y)[1]

  z <- jglog(y, alpha, lambda)
  s <- 0
  mean_spec <- apply(z, 1, mean, na.rm=TRUE)

  s <- sum((z - mean_spec)^2, na.rm=TRUE)
  
  return(s)
}

#' Performs glog transformation on the data set, using QC samples to estimate lamda and alpha parameters.
#'
#' https://doi.org/10.1186/1471-2105-8-234
#'
#' @param df Peak intensity matrix
#' @param classes Vector of class labels
#' @param qc_label Class label for QC sample
#' @export

glog_transformation <- function(df, classes, qc_label){

  if (length(which(classes %in% qc_label))==0){
    stop ("QC sample label is not present in sample class label. 
          Please check your data and qc_label parameter.")
  }
  
  df <- check_peak_matrix_orientation(peak_data = df, classes = classes)
  
  df_qc <- df[, classes == qc_label]

  # y0 gets hardocded always to 0, and later is used for alpha in glog itself
  y0 <- 0

  scal_fact <- 1
  pow_fact <- 1
  offset <- min(df_qc, na.rm=TRUE)
  df_qc <- df_qc-offset

  step_threshold <- 1e-16

  if(min(apply(df_qc, 1, var, na.rm=TRUE))==0)
  {
    varbs <- apply(df_qc, 1, var, na.rm=TRUE)
    newminVar <- sort(unique(varbs))[2]
  }else{
    newminVar <- min(apply(df_qc, 1, var, na.rm=TRUE))
  }

  low_lim <- 0
  upper_lim <- max(pmax(apply(df_qc, 1, var, na.rm=TRUE), max(apply(df_qc, 1, var, na.rm=TRUE))/newminVar))

  # What does y0 which is set to 0 is doing here?
  lambda <- optimize(f=SSE, interval=c(low_lim, upper_lim), y0, df_qc, tol=step_threshold)

  lambda <- as.numeric(lambda[[1]])

  # If SSE optimisation fails use fixed lambda value
  lambda_std <- 5.0278*10^(-09)

  error_flag <- F

  if(abs(upper_lim-lambda)<=1e-5)
  {
    cat("Error!Lambda tending to infinity!Using standard\n")
    error_flag <- TRUE
  } else if(abs(low_lim-lambda)<=1e-5)
  {
    cat("Error!Lambda tending to -infinity!Using standard\n")
    error_flag <- TRUE
  }

  if(error_flag)
  {
    lambda <- lambda_std
    scal_fact <- apply(df, 2, sum, na.rm=TRUE)
    scal_fact <- mean(scal_fact)
    scal_fact <- 1/scal_fact
  }

  df <- df*scal_fact
  df <- df^pow_fact
  df <- df-min(df, na.rm=TRUE)
  df_glog <- as.data.frame(glog(df, y0, lambda))

  return(df_glog)
}
