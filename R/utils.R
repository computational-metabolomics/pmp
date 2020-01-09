#' Return list of function argument values
#' @return List of funcation argument values in the current function call
#' @noRd

return_function_args <- function()
{
    function_argument_names <- names(formals(sys.function(sys.parent(n=1))))
    # Don't list argument values for large objects. E.g. peak matrix, class
    # labels etc.
    args_black_list <- c("df", "classes")
    function_argument_names <- function_argument_names[!function_argument_names
        %in% args_black_list]
    return(mget(function_argument_names, envir = parent.frame()))
}
