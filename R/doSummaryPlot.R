#' @importFrom gridExtra grid.arrange
#' @importFrom ggplot2 ggsave
#' 
NULL

#' Wrapper function to generate PCA plot and RSD statistics plot
#'
#' @param Data Data frame.
#' @param classes Vector of class labels.
#' @param blank Label used for blank samples, if set to NULL no samples will be removed
#' @param PQN Can be set to T or F, to perform PQN normalisation
#' @param mv_impute T or F, indicates if missing value imputation has to be carried
#' @param glogScaling T or F, applie glog transformation to the given data
#' @param scale Perform UV scaling on data
#' @param qc_label Label used for QC samples. If set to NULL, assumes that no QC samples are present in data set
#' @param ignorelabel Label for samples which should be excluded from processed data
#' @param output File name for pdf output
#' @param labels Can be set to "QC" to label only QC samples. "none" to no include labels. If set to any other value will use column names of Data.
#' @param qc_shape Shape symbol to use for QC samples
#' @param base_size Ggplot font size
#' @param pccomp PCA components to plot
#' @param plot If set to T will output pdf file, otherwise ggplot object
#' @param plotTitle Title to use for plot output
#' @return Summary plot as ggplot object or directly to the pdf file.
#' @export


doSummaryPlot <- function (Data, classes, plotTitle="PCA", blank="BLANK", PQN=F, mv_impute=T, glogScaling=T, scale=T, qc_label="QC", ignorelabel="Removed", output="PCA_plot.pdf", labels="QC", qc_shape=17, base_size = 10, pccomp=c(1,2), plot=T)
{
  Nbatches <- NULL

  if (!is.data.frame(Data) & is.list(Data))
  {
    Nbatches <- length (Data)

  }

  if (is.data.frame(Data))
  {
    Nbatches <- 1
    Data <- list(Data)
    classes <- list(classes)
    plotTitle <- list(plotTitle)
  }

  if (is.null(Nbatches))
  {
    cat ("Input data should be data frame object or list of data frames if you want to process more than one batch","\n")
    stop()
  }

  if (length (Data)!=length(classes))
  {
    cat ("Number of samples in data table doesn't match number of samples defined in \"class\" variable.")
    stop()
  }


  # RSD plot is always plotted
  plots <- vector("list", Nbatches*2)
  ind <- seq(1,length(plots),2)

  for (nb in 1:Nbatches)
  {

    sData <- data.frame(Data[[nb]])
    sclass <- classes[[nb]]
    prepData <- prepareData(Data=sData, classes=sclass, blank = blank, qc_label = qc_label, PQN=PQN, mv_impute=mv_impute, glogScaling=glogScaling, ignorelabel="Removed")

    hits4 <- which(names(prepData$RSD)==qc_label)

    if (length (hits4)>0)
    {
      QC_RSD <- prepData$RSD[[hits4]]
      QC_RSD_sum <- summary(QC_RSD)
      subt <- paste("QC RSD%:", paste(paste (names(QC_RSD_sum),round(QC_RSD_sum,1), sep=" "),collapse=", "))
      subt2 <- paste0 (" QC RSD is <30% in ", round(length(which(QC_RSD<=30))/length(QC_RSD)*100,0),"% (", length(which(QC_RSD<=30)),"/",length(QC_RSD),") of features")

    } else
    {
      subt <- NULL
      subt2 <- NULL
    }

    ## PCA
    plots[[ind[nb]]] <- doPCA (Data=prepData$Data, classes=prepData$classes, plotTitle=plotTitle[[nb]],
                               scale=scale, labels=labels, qc_shape=qc_shape, base_size = base_size,
                               pccomp=pccomp, subtitle=subt, qc_label = qc_label, PQN=PQN, mv_impute = mv_impute, glogScaling = glogScaling)

      # RSD per class
    plots[[ind[nb]+1]] <- doRSDplot (RSD=prepData$RSD, plotTitle="", base_size = 12, subtitle=subt2)

  }

  if (plot==T)
  {
    plots <- gridExtra::grid.arrange(grobs=plots, ncol=2)
    ggplot2::ggsave(filename = output, plot=plots, width=180, height = 90*Nbatches, units = "mm")
    #pdf (output, width=14, height=7*Nbatches)
    #  multiplot(plotlist=plots, layout=matrix (1:length(plots), ncol=2, byrow = T))
    #dev.off ()
    
  } else
  {
    return(plots)
  }
}

