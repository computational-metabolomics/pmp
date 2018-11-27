#' Wrapper function to transform data for statistical analysis
#'
#' @param Data Data frame.
#' @param classes Vector of class labels.
#' @param blank Label used for blank samples, if set to NULL no samples will be removed
#' @param PQN Can be set to T or F, to perform PQN normalisation
#' @param mv_impute T or F, indicates if missing value imputation has to be carried
#' @param glogScaling T or F, applie glog transformation to the given data
#' @param qc_label Label used for QC samples. If set to NULL, assumes that no QC samples are present in data set
#' @param ignorelabel Label for samples which should be excluded from processed data
#' @param checkNA removes rows or columns containing all NA's, also will chek if only all QC or analytical samples contain missing values
#' @return List of processed data table and RSD% per sample class
#' @export


prepareData <- function (Data, classes, blank="BLANK", PQN=F, mv_impute=T, glogScaling=T, qc_label="QC", ignorelabel="Removed", checkNA=T){
  
    Data <- check_peak_matrix_orientation(peak_data = Data, classes = classes)
    
    rem_samples <- NULL
    if (!is.null(blank)){
      
      rem_samples <- which(classes==blank)
    }
      
    if (!is.null(ignorelabel)){
        
      rem_samples <- append(rem_samples,which(classes==ignorelabel))
    }

    if(length(rem_samples)>0){
      
      Data <- Data[,-c(rem_samples)]
      classes <- classes[-c(rem_samples)]
    }

    # Remove rows (features) with all NA's in QC sample or in analytical sample. 
    # It can happen if we process data from more than one batch.
    
    if (checkNA==T){
      
      AllNa <- function (x) all(is.na(x))

      if (is.null(qc_label)){
        
        hits2 <- which(apply(Data,1,AllNa)==T)
      } else {
        
        hits2 <- which(apply(Data[,classes==qc_label],1,AllNa)==T)
        hits2 <- append(hits2, which(apply(Data[,-c(which(classes==qc_label))],1,AllNa)==T))
        hits2 <- unique(hits2)
      }

      if (length(hits2)>0){
        
        Data <- Data[-c(hits2),]
      }
    }
    
    RSD <- doRSD (Data=Data, classes = classes)

    #pqn normalisation
    if (PQN==T & !is.null(qc_label)){
      Data <- pqn_normalisation(df=Data, classes=classes, qc_label = qc_label)[[1]]
    }

    # mv imputation
    if (mv_impute==T){
      Data <- mv_imputation(Data, 'knn', k=5, rowmax=1, colmax=1)
    }

    #glog scaling
    if (glogScaling==T & !is.null(qc_label)){
      
      Data <- glog_transformation(df=Data, classes=classes, qc_label=qc_label)
    }

    out <- list(Data, classes, RSD)
    names (out) <- c("Data","classes","RSD")
    out
}





