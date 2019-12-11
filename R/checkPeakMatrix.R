#' @import SummarizedExperiment
#' @importFrom methods as
#' @importFrom S4Vectors DataFrame
#' @importFrom S4Vectors metadata
NULL

#' Check if peak matrix is in format features in rows, samples in columns 
#' and that all cells contain numeric values. 
#' 
#' All functions in pmp pacakge expect input peak matrix to have samples 
#' as columns and measured features in rows. This function will check input 
#' matrix orientation and will transpose it if needed. If class labels are 
#' provided this function will check if the length of labels matches number 
#' of samples in peak matrix.
#'
#' @param peak_data peak matrix
#' @param classes vector of class labels
#' @return matrix where samples are represented in columns and features in rows
#' 

check_peak_matrix <- function(peak_data, classes=NULL) {
    dims <- dim(peak_data)
    is_data_frame <- is.data.frame(peak_data)
    
    if (dims[1] < dims[2] & is.null(classes)) {
        peak_data <- t(peak_data)
        warning("Peak table was transposed to have features as rows and samples
    in columns. \n
    As there were no class labels availiable please check that peak table is \n
    still properly rotated, samples as columns and features in rows. \n
    Use 'check_df=FALSE' to keep original peak matrix orientation.")
    }
    
    if (!all(apply (peak_data, 2, typeof)!="factor" & 
        apply(peak_data, 2, typeof)!="character")){
        stop ("Peak matrix contains non-numeric values. Check your inputs!")
    }
    
    if (!is.null(classes)) {
        hits <- which(dims == length(classes))
        if (length(hits) == 2) {
            warning(" Number of samples and features is the same in your
            data matrix, 
            please make sure that you samples are in columns. \n")
        } else if (length(hits) == 0) {
            stop(" Length of sample classes doesn't match any dimension
            of input data. 
            Sample labels should match number of samples. \n")
        } else if (hits == 1) {
            # If samples are in rows, transpose data matrix
            peak_data <- t(peak_data)
        }
    }
    
    if (is_data_frame) {
        peak_data <- as.data.frame(peak_data)
    }
    
    peak_data
}

#' Check if input data is object of 'SummarizedExperiment',
#' if not convert inputs into 'SummarizedExperiment' container.
#' 
#'
#' @param peak_data peak matrix
#' @param classes vector of class labels
#' @return object of class of 'SummarizedExperiment' 
#' 
check_input_data <- function (peak_data, classes=NULL){
    if(is.null(attr(class(peak_data), "package"))){
        meta_data <- list(original_data_structure=class(peak_data))
        col_names <- colnames(peak_data)
        if (meta_data != "matrix"){
            peak_data <- as.matrix(peak_data)
        }
        peak_data <- check_peak_matrix(peak_data=peak_data,
            classes=classes)
        peak_data <- SummarizedExperiment(assays=list(peak_data))
        metadata(peak_data) <- meta_data
        if (!is.null(classes)){
            colData(peak_data) <- DataFrame(Sample=col_names, classes=classes)
        }
        colnames(peak_data) <- col_names
    }
    return(peak_data)
}

#' If input data were not of class of 'SummarizedExperiment',
#' convert output to the original R data structure type. 
#'
#' @param peak_data peak matrix
#' @return peak matrix object of the same data structure as original input data 
#' 
return_original_data_structure <- function(peak_data){
    meta_data <- metadata(peak_data)
    peak_data <- assay(peak_data)
    peak_data <- as(peak_data, meta_data$original_data_structure)
    return (peak_data)
}
