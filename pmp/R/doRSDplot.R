#' Plot violin plots from the doRSD fucntion output
#'
#' @param RSD OUtput of doRSD function
#' @param plotTitle Main tilte for the plot
#' @param base_size Font size for plot fonts
#' @param subtitle Subtitle to include in PCA plot
#' @return Ggplot object with plot(s)
#' @export

doRSDplot <- function (RSD, plotTitle=NULL, base_size = 12, subtitle=NULL)
{
  RSDs <- unlist(RSD)

  RSDsLab <- NULL
  for (rd in 1:length(RSD))
  {
    RSDsLab <- append(RSDsLab, rep(names(RSD)[rd],length(RSD[[rd]])))
  }

  plotClass <- createClassAndColors(class=RSDsLab)

  B <- data.frame (RSD=RSDs, class=plotClass$class)
  #B$class <- factor(B$class,labels = names(RSD),ordered=T)
  out <- ggplot (B, aes(y=RSD, x=class, colour=class, fill=class))+
    geom_violin(alpha=0.1, na.rm=T, show.legend = F, draw_quantiles = c(0.25,0.5,0.75))+
    scale_colour_manual(values=plotClass$manual_colors)+ theme_Publication(base_size = base_size)+ xlab("")+
    scale_y_continuous(breaks = seq(0,200,20), limits=c(0,200))+
    ggtitle(plotTitle,subtitle =  subtitle)+ ylab("RSD %")
   out
}



