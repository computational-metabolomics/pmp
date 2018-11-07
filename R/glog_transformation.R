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

#' Internal glog function to estimate SSE error
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

  df <- check_peak_matrix_orientation(peak_data = df, classes = classes)
  
  x <- df[, classes == qc_label]

  # y0 gets hardocded always to 0, and later is used for alpha in glog itself
  y0 <- 0
  N <- ncol(x)
  
  # Why max dimension here and not number of features?
  L <- max(dim(x))

  scal_fact <- 1
  pow_fact <- 1
  offset <- min(x, na.rm=TRUE)
  x <- x-offset

  step_threshold <- 1e-16

  # Why this is here if x-min(x) always will give at least one 0?
  small = min(x, na.rm=TRUE) 

  if(min(apply(x, 1, var, na.rm=TRUE))==0)
  {
    varbs <- apply(x, 1, var, na.rm=TRUE)
    newminVar <- sort(unique(varbs))[2]
  }else{
    newminVar <- min(apply(x, 1, var, na.rm=TRUE))
  }

  low_lim <- -small^2
  upper_lim <- max(pmax(apply(x, 1, var, na.rm=TRUE), max(apply(x, 1, var, na.rm=TRUE))/newminVar))

  lambda <- optimize(SSE, interval=c(low_lim, upper_lim), y0, x, tol=step_threshold)

  lambda <- as.numeric(lambda[[1]])

  # What is this?
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

  x <- df
  N <- dim(x)[2]

  if(error_flag)
  {
    lambda <- lambda_std
    scal_fact <- apply(x, 2, sum, na.rm=TRUE)
    scal_fact <- mean(scal_fact)
    scal_fact <- 1/scal_fact
  }

  x <- x*scal_fact
  x <- x^pow_fact
  x <- x-min(x, na.rm=TRUE)
  df_glog <- as.data.frame(glog(x, y0, lambda))

  return(df_glog)
}
