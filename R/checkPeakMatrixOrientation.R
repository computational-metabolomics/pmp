#' Check is input/output data are in format features in rows, samples in columns. Check if class labels are provided for all samples.
#'
#' @param df peak matrix
#' @param classes Vector of class labels
#' @return matrix in where samples are represented in columns and features in rows
#' @export

checkPeakMatrixOrientation <- function (df, classes)
{
  df <- as.matrix(df)
  dims <- dim(df)
  hits <- which(dims==length(classes))
  if (length(hits)==2)
  {
    warning (" Number of samples and features is the same in your data matrix, please make sure that you samples are in columns. \n")
  }
  
  if (length(hits)==0)
  {
    stop (" Length of sample classes doesn't match any dimension of input data. Sample labels should match number of samples. \n")
  }
  
  # If samples are in rows, transpose data matrix
  if (hits==1)
  {
    df <- t(df)
  }
  
  df
}