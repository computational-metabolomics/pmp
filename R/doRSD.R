#' Calculate RSD\% values per sample group
#'
#' @param Data Data frame
#' @param classes Vector of class labels
#' @return List of RSD\% values for each variable for each sample group
#' @export

doRSD <- function(Data, classes){
  
  Data <- check_peak_matrix_orientation(peak_data = Data, classes = classes)
  
  cl <- unique (classes)
  out <- vector("list",length(cl))
  names (out) <- cl

  #Calculate RSD before scaling and MV imputation
  FUN = function(x) sd(x,na.rm=T)/mean(x,na.rm=T)*100.0

  for (slab in 1: length(out))
  {
    out[[slab]] <- apply(Data[,classes==names(out)[slab]],1,FUN)
  }
  out
}
