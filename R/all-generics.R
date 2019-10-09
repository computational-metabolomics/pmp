setGeneric(name="setScaledPeakMatrix", 
    def=function(glog_object, scaled_peak_matrix){
        standardGeneric("setScaledPeakMatrix") 
})

setGeneric(name="setGlogLambdaSummary", 
    def=function(glog_object, lambda, lambda_opt, error_flag){
        standardGeneric("setGlogLambdaSummary") 
})

setGeneric(name="setGlogLambdaOptimisationPlot", 
    def=function(glog_object, lambda_optimisation_plot){
        standardGeneric("setGlogLambdaOptimisationPlot") 
})

setGeneric(name="getScaledPeakMatrix", 
    def=function(glog_object){
        standardGeneric("getScaledPeakMatrix") 
})

setGeneric(name="getGlogLambdaOptimisationPlot", 
    def=function(glog_object){
        standardGeneric("getGlogLambdaOptimisationPlot") 
})