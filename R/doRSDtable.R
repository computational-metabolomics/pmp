
#' Plot violin plots from the doRSD fucntion output
#'
#' @param RSD OUtput of doRSD function
#' @param QC_label Label used for QC samples. If set to NULL, assumes that no QC samples are present in data set
#' @param Blank_label Label used for Blank samples
#' @return Table of RSD% values per group
#' @export

doRSDtable <- function (RSD, QC_label="QC", Blank_label="Blank")
{
  tableCNames <- names(RSD)

  #tableCNames <- reorderFactorLevels(labels = tableCNames, QC_label = QC_label, Blank_label = Blank_label)
  tableCNames <- createClassAndColors(class = tableCNames, QC_label = QC_label, Blank_label = Blank_label)$class

  tableData <- lapply(RSD, summary)

  # So summary table doesn't return NA count if there are no NA's. Honestly?
  for (i in 1:length(tableData))
  {
    if (length(tableData[[i]])<7)
    {
      tableData[[i]] <- append(tableData[[i]],0)
    }
  }

  tableData <- do.call(rbind, tableData)
  colnames(tableData)[7] <- "NA's"

  tableData <- as.data.frame(tableData)

  rownames(tableData) <- names (RSD)

  tableData <- tableData[order(factor(rownames(tableData),levels=tableCNames,ordered=T)),]
  tableData
}



