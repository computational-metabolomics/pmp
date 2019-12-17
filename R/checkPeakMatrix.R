#' @import SummarizedExperiment
#' @importFrom methods as
#' @importFrom methods is
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
    
    if (dims[1] < dims[2] & is.null(classes)) {
        peak_data <- t(peak_data)
        warning("Peak table was transposed to have features as rows and samples
    in columns. \n
    As there were no class labels availiable please check that peak table is \n
    still properly rotated, samples as columns and features in rows. \n
    Use 'check_df=FALSE' to keep original peak matrix orientation.")
    }
    
    if (!is.numeric(as.matrix(peak_data))){
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
    meta_data <- list(original_data_structure=class(peak_data)[1])
    
    if(is(peak_data, "SummarizedExperiment")){
        metadata(peak_data)$original_data_structure <- "SummarizedExperiment"
    
    } else {

        if (meta_data$original_data_structure != "matrix"){
            peak_data <- as.matrix(peak_data)
        }
        peak_data <- check_peak_matrix(peak_data=peak_data,
            classes=classes)
        col_names <- colnames(peak_data)
        peak_data <- SummarizedExperiment(assays=list(peak_data))
        metadata(peak_data) <- meta_data
        if (!is.null(classes)){
            colData(peak_data) <- DataFrame(classes=classes)
        }
        colnames(peak_data) <- col_names
    }
    return(peak_data)
}

#' If input data were not of class of 'SummarizedExperiment',
#' convert output to the original R data structure type. All values stored in 
#' 'metadata' slot of \link[SummarizedExperiment]{SummarizedExperiment} object 
#' will be exported as 'attributes' of the output object.
#'
#' @param summarized_experiment_object peak matrix
#' @return peak matrix object of the same data structure as original input data 
#' 
return_original_data_structure <- function(summarized_experiment_object){
    meta_data <- metadata(summarized_experiment_object)
    if (!meta_data$original_data_structure == "SummarizedExperiment"){
        peak_data <- assay(summarized_experiment_object)
        # as() can't convert matrix to data.frame, but works with all other objects
        if (meta_data$original_data_structure == "data.frame"){
            peak_data <- as.data.frame(peak_data)
        } else if (meta_data$original_data_structure != "matrix"){
            peak_data <- as(peak_data, meta_data$original_data_structure)
        }
    
        ## Add all metadata as output object attributes, not sure if this works for
        ## DataFrame class as well.
        meta_data$original_data_structure <- NULL
        attributes(peak_data) <- c(attributes(peak_data), meta_data)
        return (peak_data)
    } else {
        return (summarized_experiment_object)
    }
}
