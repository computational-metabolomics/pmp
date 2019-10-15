setGeneric(name="setGlogScaledPeakMatrix", 
    def=function(glog_object, scaled_peak_matrix){
        standardGeneric("setGlogScaledPeakMatrix") 
})

setGeneric(name="setGlogLambdaSummary", 
    def=function(glog_object, lambda, lambda_opt, error_flag){
        standardGeneric("setGlogLambdaSummary") 
})

setGeneric(name="setGlogLambdaOptimisationPlot", 
    def=function(glog_object, lambda_optimisation_plot){
        standardGeneric("setGlogLambdaOptimisationPlot") 
})

setGeneric(name="getGlogScaledPeakMatrix", 
    def=function(glog_object){
        standardGeneric("getGlogScaledPeakMatrix") 
})

setGeneric(name="getGlogLambdaOptimisationPlot", 
    def=function(glog_object){
        standardGeneric("getGlogLambdaOptimisationPlot") 
})