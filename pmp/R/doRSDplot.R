doRSDplot <- function (RSD, plotTitle="PCA", base_size = 12, subtitle=NULL)
{
  RSDs <- unlist(RSD)

  RSDsLab <- NULL
  for (rd in 1:length(RSD))
  {
    RSDsLab <- append(RSDsLab, rep(names(RSD)[rd],length(RSD[[rd]])))
  }

  B <- data.frame (RSD=RSDs, class=RSDsLab)
  out <- ggplot (B, aes(y=RSD, x=class, colour=class, fill=class))+
    geom_violin(alpha=0.1, na.rm=T, show.legend = F)+
    scale_colour_Publication()+ theme_Publication(base_size = base_size)+ xlab("")+
    scale_y_continuous(breaks = seq(0,200,20), limits=c(0,200))+
    ggtitle("",subtitle =  subtitle)+ ylab("RSD %")
   out
}



