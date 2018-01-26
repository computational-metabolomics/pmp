doPCA <- function (Data, classes, plotTitle="PCA",PQN, mv_impute, glogScaling, scale=F, labels="QC", qc_label, qc_shape=17, base_size = 12, pccomp=c(1,2), subtitle=NULL)
{
  require (ggplot2)
  require (reshape2)
  require (gridExtra)
  require (ggthemes)
  # ggplot publication theme.
  source ("//its-rds/2015/viantm-01/users/jankevia/Github/pmp/pmp/R/ggplot_theme_pub.R")
  source ("//its-rds/2015/viantm-01/users/jankevia/Github/pmp/pmp/R/multiplot.R")

  ## PCA
  PCA <- prcomp(t(Data), center=T, scale=scale)
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

  A <- data.frame (class=as.factor(classes),x=PCA$x[,pccomp[1]], y=PCA$x[,pccomp[2]], labels=slabels, shapes=shapes)
  labx <- paste("PC",pccomp[1]," ",varPCA[pccomp[1]],"%", sep="")
  laby <- paste("PC",pccomp[2]," ",varPCA[pccomp[2]],"%", sep="")

  labx <- paste0(labx," (",paste("PQN:",PQN,", glog scaling:",glogScaling, ", UV scaling:",scale,sep=""),")")

   out <- ggplot (data=A, aes(x=x, y=y, color=class, label=labels, shape=shapes))+
      geom_point(na.rm=T)+ scale_shape_identity()+
      xlab (labx) + ylab (laby)+
      ggtitle(plotTitle, subtitle=subtitle)+
      geom_text(na.rm=T)+
      stat_ellipse()+
      scale_colour_Publication()+ theme_Publication(base_size = base_size)

      # Looks like a bug in ggplot, label always appear in bold and weirdly scaled. Doesn't respond to fontface and family parameters
      #geom_text (x=Inf, y=-Inf, label=annotateText, hjust=1, vjust=-0.5, colour="black", size=3, fontface="plain", family="Times", lineheight=0.8)
      #annotate ("text", x=Inf, y=-Inf, label=annotateText, size=3, colour="black")

   out
}



