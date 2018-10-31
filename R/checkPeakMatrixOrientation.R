#' Check is input/output data are in format features in rows, samples in columns. Check if class labels are provided for all samples.
#'
#' @param peak_data peak matrix
#' @param classes Vector of class labels
#' @return matrix where samples are represented in columns and features in rows
#' @export

check_peak_matrix_orientation <- function (peak_data, classes=NULL)
{
  dims <- dim(peak_data)
  if (dims[1] < dims[2] & is.null(classes)){
    peak_data <- t(peak_data)
    warning ("Peak table transposed to have features as rows and samples in columns. \n
             As there were no class labels availiable please check that you peak table is
             still properly rotated.")
  }
  
  if (!is.null(classes)){
    hits <- which(dims==length(classes))
    if (length(hits)==2){
      warning (" Number of samples and features is the same in your data matrix, 
               please make sure that you samples are in columns. \n")
    }
  
    if (length(hits)==0){
      stop (" Length of sample classes doesn't match any dimension of input data. 
            Sample labels should match number of samples. \n")
    }
  
    # If samples are in rows, transpose data matrix
    if (hits==1){
      peak_data <- t(peak_data)
    }
  }
  
  peak_data
}
