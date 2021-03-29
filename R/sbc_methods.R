#' @importFrom stats smooth.spline predict
#' @importFrom utils capture.output
NULL

##' Transform data using natural logarithm
##'
##' @param x Values to convert.
##' @param a log offset.
##' @param inverse To apply log or reverse back to original values.
##' 
##' @return vector of log transformed values
##' @noRd

logFit <- function (x, a=1, inverse=FALSE) {
    if (inverse) {
        out <- 0.25 * exp(-x) * (4 * exp(2 * x) - (a * a))
    }else{
        out <- log((x + sqrt(x^2 + a^2))/2)
    }
    return(out)
}

##' Cubic smoothing spline fitting function for batch effect correction
##'
##' @param x QC sample measurement order indices.
##' @param y Signal intensities.
##' @param newX Range of sample indices to predict.
##' @param log TRUE or FALSE to perform the signal correction fit on the log
##'scaled data. Default is TRUE.
##' @param a log scaling offset.
##' @param spar Smoothing parameter of the fitted curve. Should be in the range
##'0 to 1. If set to 0 it will be estimated using leave-one-out
##'cross-validation.
##' @return vector of cubic smoothing spline fitted values
##' @noRd

splineSmoother <- function(x, y, newX, log=log, a=1, spar) {
    if(log == TRUE) {
        y <- logFit(y, a=a)
    }
    if (spar == 0) {
        # Supress spline fitting CV messages cluttering terminal output
        smoothSplineMessages <- capture.output(sp.obj 
            <- smooth.spline(x, y, cv=TRUE), file=NULL, type="message")
    } else {
        smoothSplineMessages <- character()
        sp.obj <- smooth.spline(x, y, spar=spar)
    }

    valuePredict <- predict(sp.obj, newX)

    if(log==TRUE) {
        valuePredict$y <- logFit(valuePredict$y, a=a, inverse=TRUE)
    }
    
    return(valuePredict$y)
}


##' Signal batch correction wrapper function
##'
##' @param id Feature id.
##' @param qcData Data frame with QC samples.
##' @param qcBatch QC batch numbers.
##' @param qcOrder QC measurement order.
##' @param spar Spline smoothing parameter. Should be in the range 0 to 1.
##'If set to 0 it will be estimated using leave-one-out cross-validation.
##' @param log TRUE or FALSE to perform the signal correction fit on the
##'log scaled data. Default is TRUE.
##' @param batch A vector indicating the batch each sample was measured in.
##'If only one batch was measured then all values should be set to 1.
##' @param minQC Minimum number of QC samples required for signal correction.
##' @param order A numeric vector indicating the order in which samples
##'were measured.
##' @return vector of corrected values of selected feature for QC data
##' @noRd

sbcWrapper <- function(id, qcData, order, qcBatch, qcOrder, log=log, spar=spar,
    batch=batch, minQC) {
    
    out <- tryCatch ({

    # qc data for feature
    subData <- qcData[id, ]
    
    # unique batches
    nbatch <- unique(qcBatch)

    # preallocate output
    outl <- numeric(length(order))

    # for each batch
    for (nb in seq_len(length(nbatch))) {
        # get the injection order for the QCs
        x <- qcOrder[qcBatch == nbatch[nb]]
        # get the response for the QCs
        y <- subData[qcBatch == nbatch[nb]]
        
        # remove any index with NA in the response
        NAhits <- which(is.na(y))

        if (length(NAhits) > 0) {
            x <- x[-c(NAhits)]
            y <- y[-c(NAhits)]
        }

        # if we have enough QCs values
        if (length(y) >= minQC) {
            # fit spline to QCs and get predictions for all samples in batch
            outl[batch==nbatch[nb]] <- splineSmoother(x=x, y=y, newX=order[batch==nbatch[nb]], log=log,
            a=1, spar=spar)
        } else {
            # otherwise replace with NA
            outl[batch==nbatch[nb]] = NA
        }
    }

    outp=outl
    },

    error = function(e){
        print (e, "\n")
        return(NULL)
    })
    return(outp)
}
