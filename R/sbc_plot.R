#' @import ggplot2
#' @importFrom grDevices pdf dev.off
#' @importFrom reshape2 melt
NULL

#' Plot QCRSC corrected outputs
#' 
#' Plot the output from signal batch correction for the selected or
#' the first 100 features.
#'
#' @inheritParams filter_peaks_by_blank
#' @param corrected_df Output from \link[pmp]{QCRSC} function.
#' @param batch \code{numeric()} or \code{character()}, a vector indicating
#' the batch each sample was measured in. If only one batch was measured then
#' all values should be set to 1
#' @param indexes \code{numeric()}, a vector defining which features to plot.
#' If set to \code{NULL} will plot the first 100.
#' @param output \code{character()}, a filename of the output pdf file. 
#' Can include the path. If set to \code{NULL} output will be list object 
#' containing class \code{ggplot} plots.
#' 
#' @return Pdf file or \code{list()} object \code{ggplot} class showing data 
#' before and after signal correction.
#' 
#' @examples 
#' 
#' order <- c(1:ncol(MTBLS79))
#' data <- MTBLS79[1:10, ]
#' 
#' out <- QCRSC(df =data, order=order, batch=MTBLS79$Batch, 
#'    classes=MTBLS79$Class, spar=0, minQC=4)
#' plots <- sbc_plot (df=data, corrected_df=out, classes=MTBLS79$Class,
#'    batch=MTBLS79$Batch, output=NULL) 
#'
#' @export

sbc_plot <- function(df, corrected_df, classes, batch, indexes = NULL,
    qc_label="QC", output = "sbcms_plots.pdf") {
    df <- check_input_data(df=df, classes=classes)
    corrected_df <- check_input_data(df=corrected_df, classes=classes)
    shapes <- rep(19, length(classes))
    shapes[classes == qc_label] <- 3
    manual_color <- c("#386cb0", "#ef3b2c", "#7fc97f", "#fdb462", "#984ea3", 
        "#a6cee3", "#778899", "#fb9a99", "#ffff33")
    gg_THEME <- theme(panel.background = element_blank(), 
        panel.grid.major = element_line(color = "gray80", linewidth = 0.3), 
        axis.line = element_line(color = "black"), 
        axis.text = element_text(color = "black"), 
        axis.title = element_text(color = "black"), 
        panel.grid.minor.x = element_line(color = "gray80", 
            linewidth = 0.3, linetype = "dashed"), 
        panel.grid.minor.y = element_line(color = "gray80", linewidth = 0.3))
    plots <- list()
    if (is.null(indexes) & nrow(df) >= 100) {
        indexes <- seq_len(100)
    } else if (nrow(df) < 100 & is.null(indexes)) {
        indexes <- seq_len(nrow(df))
    }
    for (peakn in indexes) {
        A <- data.frame(x = c(seq_len(ncol(df))), 
            original = as.vector(log(assay(df[peakn, ]), 10)),
            corrected = as.vector(log(assay(corrected_df[peakn, ]), 10)), 
            batch = as.factor(batch), shapes = shapes)
        A <- melt(A, id.vars = c("x", "batch", "shapes"))
        
        plots[[peakn]] <- ggplot(A, 
            aes(x = .data[['x']], 
                y = .data[['value']], 
                col = .data[['batch']], 
                shape = .data[['shapes']])) + 
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
