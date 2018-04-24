#' @import stats
#' @import ggplot2
NULL

#' Perform PCA analysis ot single data set of list of data sets
#'
#' @param Data Data frame.
#' @param classes Vector of class labels.
#' @param plotTitle Character value to display as the plot title
#' @param PQN Can be set to T or F, to perform PQN normalisation
#' @param mv_impute T or F, indicates if missing value imputation has to be carried
#' @param glogScaling T or F, applie glog transformation to the given data
#' @param scale Perform UV scaling on data
#' @param labels Can be set to "QC" to label only QC samples. "none" to no include labels. If set to any other value will use column names of Data.
#' @param qc_label Label used for QC samples. If set to NULL, assumes that no QC samples are present in data set
#' @param qc_shape Shape symbol to use for QC samples
#' @param base_size Font size for plot fonts
#' @param pccomp PCA components to plot
#' @param subtitle Subtitle to include in PCA plot
#' @param loadings T or F, to include PCA loadings plot or not
#' @param loadingsCol Colors to use for loadings plot
#' @return Ggplot object with plot(s)
#' @export

doPCA <- function (Data, classes, plotTitle="PCA",PQN, mv_impute, glogScaling, scale=F, labels="QC", qc_label, qc_shape=17, base_size = 12, pccomp=c(1,2), subtitle=NULL,
                   loadings=F, loadingsCol=NULL)
{
  ## PCA
  PCA <- stats::prcomp(t(Data), center=T, scale=scale)
  varPCA <- round(PCA$sdev^2/sum(PCA$sdev^2)*100,1)

  slabels=colnames(Data)
  # Make QC samples to be labeled
  shapes <- rep(19,length(classes))
  shapes[classes==qc_label] <- qc_shape

  # Labels
  if (labels=="QC")
  {
    sh <- which(classes==qc_label)
    slabels[-c(sh)] <- ""
    rm (sh)
  }

  if (labels=="none")
  {
      slabels <- rep(NA,length(slabels))
  }

  plotClass <- createClassAndColors(class=classes)


  A <- data.frame (class=plotClass$class, pc1=PCA$x[,pccomp[1]], pc2=PCA$x[,pccomp[2]], labels=slabels, shapes=shapes)
  labx <- paste("PC",pccomp[1]," ",varPCA[pccomp[1]],"%", sep="")
  laby <- paste("PC",pccomp[2]," ",varPCA[pccomp[2]],"%", sep="")

  labx <- paste0(labx," (",paste("PQN:",PQN,", glog scaling:",glogScaling, ", UV scaling:",scale,sep=""),")")

   out <- ggplot2::ggplot (data=A, aes(x=pc1, y=pc2, color=class, label=labels, shape=shapes))+
      geom_point(na.rm=T)+ scale_shape_identity()+
      xlab (labx) + ylab (laby)+
      ggtitle(plotTitle, subtitle=subtitle)+
      geom_text(na.rm=T)+
      stat_ellipse()+
     scale_colour_manual(values=plotClass$manual_colors)+ theme_Publication(base_size = base_size)

      # Looks like a bug in ggplot, label always appear in bold and weirdly scaled. Doesn't respond to fontface and family parameters
      #geom_text (x=Inf, y=-Inf, label=annotateText, hjust=1, vjust=-0.5, colour="black", size=3, fontface="plain", family="Times", lineheight=0.8)
      #annotate ("text", x=Inf, y=-Inf, label=annotateText, size=3, colour="black")


   if (loadings==T)
   {
     if (is.null(loadingsCol)) loadingsCol <- rep(1,nrow(Data))
     B <- data.frame (pc1=PCA$rotation[,pccomp[1]], pc2=PCA$rotation[,pccomp[2]], labels=rownames(Data), loadingsCol=loadingsCol)

     out2 <- ggplot2::ggplot (data=B, aes(x=pc1, y=pc2, color=1, label=labels, shape=3))+
       scale_shape_identity()+
       xlab ("Loadings, PC1") + ylab ("Loadings, PC2")+
       geom_text(na.rm=T, hjust=0, color=loadingsCol)+
       theme_Publication(base_size = base_size)+ theme(legend.position="none")
     out <- list(out,out2)
   }
   out
}



