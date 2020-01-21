#' Return list of function argument values
#' @return List of function argument values in the current function call
#' @noRd

return_function_args <- function()
{
    function_argument_names <- names(formals(sys.function(sys.parent(n=1))))
    # Don't list argument values for large objects. E.g. peak matrix, class
    # labels etc.
    args_black_list <- c("df", "classes", "ref_mean")
    function_argument_names <- function_argument_names[!function_argument_names
        %in% args_black_list]
    return(mget(function_argument_names, envir = parent.frame()))
}

#' Return history of applied functions and argument from pmp package.
#' 
#' @inheritParams filter_peaks_by_blank
#' @return List of function  names and argument values.
#' 
#' @examples
#' df <- MTBLS79[ ,MTBLS79$Batch == 1]
#' df$Class[1:2] <- "Blank"
#' out <- filter_peaks_by_blank(df=df, fold_change=1.2, 
#'    classes=df$Class, blank_label="Blank", qc_label=NULL, 
#'    remove_samples=FALSE, remove_peaks=TRUE, fraction_in_blank=0)
#' processing_history(out)
#' 
#' @export

processing_history <- function(df){
    if(is(df, "SummarizedExperiment")){
        processing_history <- metadata(df)$processing_history
    } else {
        processing_history <- attributes(df)$processing_history
    }
        
    return (processing_history)
}