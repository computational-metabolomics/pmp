prepareData <- function (Data, classes, blank="BLANK", PQN=F, mv_impute=T, glogScaling=T, qc_label="QC", ignorelabel="Removed")
{
    shits <- NULL
    if (!is.null(blank))
    {
      shits <- which(classes==blank)
    }
      if (!is.null(ignorelabel))
    {
      shits <- append(shits,which(classes==ignorelabel))
    }

    if(length(shits)>0)
    {
      Data <- Data[,-c(shits)]
      classes <- classes[-c(shits)]
    }

    # Remove rows with all NA's in QC sample or in analytical sample, this can happen if we process data from more than one batch.
    AllNa <- function (x) all(is.na(x))

    hits2 <- which(apply(Data[,classes==qc_label],1,AllNa)==T)
    hits2 <- append(hits2, which(apply(Data[,-c(which(classes==qc_label))],1,AllNa)==T))
    hits2 <- unique(hits2)

    if (length(hits2)>0)
    {
      Data <- Data[-c(hits2),]
    }

    RSD <- doRSD (Data=Data, classes = classes)

    #pqn normalisation
    if (PQN==T)
    {
      source ("//its-rds/2015/viantm-01/users/jankevia/Github/pmp/pmp/R/normalisation.R")
      Data <- pqn_normalisation(df=Data, classes=classes, qc_label = qc_label)[[1]]
    }

    if (mv_impute==T)
    {
      #impute missing values
      source ("//its-rds/2015/viantm-01/users/jankevia/Github/pmp/pmp/R/mv_imputation.R")
      vals <- mv_imputation(t(Data), 'knn', k=5, rowmax=1, colmax=1)
      Data <- t(vals)
    }

    #glog scaling
    if (glogScaling==T)
    {
      source("//its-rds/2015/viantm-01/users/jankevia/Github/pmp/pmp/R/glog_transformation.R")
      Data <- glog_transformation(df=t(Data), classes=classes, qc_label=qc_label)
      Data <- t(Data)
    }

    out <- list(Data, classes, RSD)
    names (out) <- c("Data","classes","RSD")
    out
}



