#' @importFrom methods new
#'
NULL

#' S4 class for Variance stabilising generalised lograrithm transformation
#'  
#' This class defines structure and methods of glog_transformation output.
#' 
#' @slot scaled_peak_matrix a data frame, peak inntensity matrix after glog
#' transformation
#' @slot lambda_summary a list containing optimised lambda value and if 
#' optimisation has failed
#' @slot lambda_optimisation_plot ggplot2 object of glog lambda optimisation
#' @export

GlogOutput <- setClass(
    "GlogOutput",
    slots=c(scaled_peak_matrix="data.frame",
        lambda_summary="list",
        lambda_optimisation_plot="list")
)

#' Method to set peak matrix slot after glog scaling
#' 
#' @aliases setScaledPeakMatrix
#' @param glog_object R object of 'GlogOuput' class
#' @param scaled_peak_matrix a data frame, peak inntensity matrix after glog
#' transformation
#' @export
#' 
setMethod(f="setScaledPeakMatrix", signature="GlogOutput",
    definition=function(glog_object, scaled_peak_matrix){
    glog_object@scaled_peak_matrix <- scaled_peak_matrix
    return (glog_object)
})

#' Method to set glog lambda optimisation summary
#' 
#' @aliases setGlogLambdaSummary
#' @param glog_object R object of 'GlogOuput' class
#' @param lambda glog lambda value
#' @param lambda_opt glog optimised lambda value
#' @param error_flag logical value from SSE optimisation call
#' transformation
#' 
setMethod(f="setGlogLambdaSummary", signature="GlogOutput",
    definition=function(glog_object, lambda, lambda_opt, error_flag){
    lambda_summary <- list (lambda=lambda, lambda_opt=lambda_opt,
        error_flag=error_flag)
    glog_object@lambda_summary <- lambda_summary
    return (glog_object)
})

#' Method to set lamba optimisation plot slot
#' 
#' @aliases setGlogLambdaOptimisationPlot
#' @param glog_object R object of 'GlogOuput' class
#' @param lambda_optimisation_plot output from 'glog_plot_optimised_labmda'
#' transformation
#' 
setMethod(f="setGlogLambdaOptimisationPlot", signature="GlogOutput",
    definition=function(glog_object, lambda_optimisation_plot){
    glog_object@lambda_optimisation_plot <- list(lambda_optimisation_plot)
    return (glog_object)
})
