#' @import ggplot2
#' @importFrom grDevices pdf dev.off
#' @importFrom reshape2 melt
NULL

##' Plot the output from signal batch correction for the selected or
##'the first 100 features
##'
##' @param df A data frame containing the original data before correction
##'(samples in columns, features in rows).
##' @param corrected_df A data frame containing the corrected data,
##'as output by QCRSC
##' @param batch A vector indicating the batch each sample was measured in.
##'If only one batch was measured then all values should be set to 1
##' @param classes A factor or character vector of sample classes.
##'All QC samples should be labelled 'QC'.
##' @param indexes Numeric vector defining whihc features from data frame to 
##'plot. If set to NULL will plot the first 100.
##' @param output Filename of the output pdf file. Can include the path.
##'If set to NULL output will be list object containing ggplot2 plots.
##' @param qc_label Class label for QC sample.
##' 
##' @return Pdf file or list object showing data before and after 
##'signal correction
##' 
##' @examples 
##' 
##' classes <- MTBLS79$Class
##' batch <- MTBLS79$Batch
##' order <- c(1:ncol(MTBLS79))
##' data <- SummarizedExperiment::assay(MTBLS79[1:10, ] )
##' 
##' out <- QCRSC(d =data, order=order, batch=batch, classes=classes,
##'     spar=0, minQC=4)
##' plots <- sbcmsPlot (df=data, corrected_df=out, classes, batch,
##'     output=NULL) 
##'
##' @export

sbcmsPlot <- function(df, corrected_df, classes, batch, indexes = NULL,
    qc_label="QC", output = "sbcms_plots.pdf") {
    
    shapes <- rep(19, length(classes))
    shapes[classes == qc_label] <- 3
    manual_color <- c("#386cb0", "#ef3b2c", "#7fc97f", "#fdb462", "#984ea3", 
        "#a6cee3", "#778899", "#fb9a99", "#ffff33")
    gg_THEME <- theme(panel.background = element_blank(), 
        panel.grid.major = element_line(color = "gray80", size = 0.3), 
        axis.line = element_line(color = "black"), 
        axis.text = element_text(color = "black"), 
        axis.title = element_text(color = "black"), 
        panel.grid.minor.x = element_line(color = "gray80", 
        size = 0.3, linetype = "dashed"), 
        panel.grid.minor.y = element_line(color = "gray80", size = 0.3))
    plots <- list()
    
    if (is.null(indexes) & nrow(df) >= 100) {
        indexes <- seq_len(100)
    } else if (nrow(df) < 100 & is.null(indexes)) {
        indexes <- seq_len(nrow(df))
    }
    for (peakn in indexes) {
        A <- data.frame(x = c(seq_len(ncol(df))), 
            original = log(df[peakn, ], 10),
            corrected = log(corrected_df[peakn, ], 10), 
            batch = as.factor(batch), shapes = shapes)
        A <- melt(A, id.vars = c("x", "batch", "shapes"))
        
        plots[[peakn]] <- ggplot(A, 
            aes_(~x, ~value, col = ~batch, shape = ~shapes)) + 
            facet_grid(variable ~ .) + geom_point() + 
            scale_shape_identity() + geom_point(size=2) +
            scale_colour_manual(values = manual_color) + 
            ggtitle(row.names(df)[peakn]) + ylab("log10(intensity)") + 
            xlab("injection order") + gg_THEME
    }
    # remove lists with NULL values
    plots <- plots[lengths(plots) != 0]
    
    if (!is.null(output)) {
        pdf(output)
        invisible(lapply(plots, print))
        dev.off()
    } else {
        return(plots)
    }
}
