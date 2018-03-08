doRSDtable <- function (RSD, QC_label="QC", Blank_label="Blank")
{
  tableCNames <- names(RSD)

  tableCNames <- reorderFactorLevels(labels = tableCNames, QC_label = QC_label, Blank_label = Blank_label)


  tableData <- lapply(RSD, summary)

  # So summary table doesn't return NA count of there are no NA's. Honestly?
  for (i in 1:length(tableData))
  {
    if (length(tableData[[i]])<7)
    {
      tableData[[i]] <- append(tableData[[i]],0)
    }
  }

  tableData <- do.call(rbind, tableData)
  colnames(tableData)[7] <- "NA's"

  rownames(tableData) <- names (RSD)

  tableData <- tableData[order(factor(rownames(tableData),levels=tableCNames,ordered=T)),]
  tableData
}



