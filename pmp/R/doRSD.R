doRSD <- function(Data, classes)
{
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
