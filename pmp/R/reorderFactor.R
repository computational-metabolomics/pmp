reorderFactorLevels <- function (labels, QC_label="QC", Blank_label="Blank")
{

  reorderNames <- sort(as.character(unique(labels)))

  hit1 <- which(reorderNames==QC_label)
  if (length(hit1)==0) hit1 <- NULL

  hit2 <- which(reorderNames==Blank_label)
  if (length(hit2)==0) hit2 <- NULL

  if (!is.null(c(hit1,hit2)))
  {
    remo <- c(1:length(reorderNames))[-c(hit1,hit2)]
  } else
  {
    remo <- c(1:length(reorderNames))
  }

  reorderNames <- reorderNames[c(hit1, hit2,remo)]
  reorderNames
}