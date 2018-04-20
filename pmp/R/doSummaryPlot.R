doSummaryPlot <- function (Data, classes, plotTitle="PCA", blank="BLANK", PQN=F, mv_impute=T, glogScaling=T, scale=T, qc_label="QC", ignorelabel="Removed", output="PCA_plot.pdf", labels="QC", qc_shape=17, base_size = 12, pccomp=c(1,2), plot=T)
{
  require (ggplot2)
  require (reshape2)
  require (gridExtra)
  require (ggthemes)
  # ggplot publication theme.

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
    pdf (output, width=14, height=7*Nbatches)
      multiplot(plotlist=plots, layout=matrix (1:length(plots), ncol=2, byrow = T))
    dev.off ()
  } else
  {
    return(plots)
  }
}

