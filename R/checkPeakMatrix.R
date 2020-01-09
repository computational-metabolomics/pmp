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
#' @param df A matrix-like (e.g. and ordinary matrix, a data frame) object with
#' all values of class \code{numeric()} or \code{integer()} of peak
#' intensities, areas or other quantitative characteristic.
#' @param classes \code{character()}, vector of class labels. Must be the same 
#' length as the number of sample in the input peak table.
#' @return \code{numeric()}, matrix-like object where samples are represented
#' in columns and features in rows.
#' 
#' @noRd

check_peak_matrix <- function(df, classes=NULL) {
    dims <- dim(df)
    if (dims[1] < dims[2] & is.null(classes)) {
        df <- t(df)
        warning("Peak table was transposed to have features as rows and samples
    in columns. \n
    As there were no class labels availiable please check that peak table is \n
    still properly rotated, samples as columns and features in rows. \n
    Use 'check_df=FALSE' to keep original peak matrix orientation.")
    }
    if (!is.numeric(as.matrix(df))){
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
            df <- t(df)
        }
    }
    df
}

#' If needed convert input data to object of class 
#' \link[SummarizedExperiment]{SummarizedExperiment}
#' 
#' @param df A matrix-like (e.g. an ordinary matrix, a data frame) object with
#' all values of class \code{numeric()} or \code{integer()} of peak
#' intensities, areas or other quantitative characteristic.
#' @param classes \code{character()}, vector of class labels. Must be the same 
#' length as the number of sample in the input peak table.
#' @return object of class of 
#' \link[SummarizedExperiment]{SummarizedExperiment}. 
#' 
#' @noRd
check_input_data <- function (df, classes=NULL){
    meta_data <- list(original_data_structure=class(df)[1])
    if(is(df, "SummarizedExperiment")){
        metadata(df)$original_data_structure <- "SummarizedExperiment"
    } else {
        if (meta_data$original_data_structure != "matrix"){
            df <- as.matrix(df)
        }
        df <- check_peak_matrix(df=df,
            classes=classes)
        col_names <- colnames(df)
        df <- SummarizedExperiment(assays=list(df))
        metadata(df) <- meta_data
        if (!is.null(classes)){
            colData(df) <- DataFrame(classes=classes)
        }
        colnames(df) <- col_names
    }
    return(df)
}

#' If input data object was not of class of 
#' \link[SummarizedExperiment]{SummarizedExperiment},
#' convert output to the original R data structure. 
#' 
#' All values stored in \code{metadata} element of \code{SummarizedExperiment}
#' will be exported as \link[base]{attributes} of the output object.
#'
#' @param summarized_experiment_object object of class
#' \code{SummarizedExperiment}
#' @return A matrix-like (e.g. an ordinary matrix, a data frame) object with
#' all values of class \code{numeric()}. 
#' Values from \code{colData} and \code{rowData} elements are not returned. 
#' 
#' @noRd
return_original_data_structure <- function(summarized_experiment_object){
    meta_data <- metadata(summarized_experiment_object)
    if (!meta_data$original_data_structure == "SummarizedExperiment"){
        peak_data <- assay(summarized_experiment_object)
        # as() can't convert matrix to data.frame, but works with all other 
        # objects
        if (meta_data$original_data_structure == "data.frame"){
            peak_data <- as.data.frame(peak_data)
        } else if (meta_data$original_data_structure != "matrix"){
            peak_data <- as(peak_data, meta_data$original_data_structure)
        }
        # Add all metadata as output object attributes
        meta_data$original_data_structure <- NULL
        attributes(peak_data) <- c(attributes(peak_data), meta_data)
        if (ncol(rowData(summarized_experiment_object)) != 0){
            peak_data <- list(df=peak_data, 
                flags=rowData(summarized_experiment_object))
        }
        return (peak_data)
    } else {
        return (summarized_experiment_object)
    }
}
