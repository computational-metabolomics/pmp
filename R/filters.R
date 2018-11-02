#' Blank filter
#'
#' @param df Peak intensity matrix
#' @param fold_change Minimum fold change between analytical and blank samples.
#' @param classes Vector of class labels
#' @param blank_label Class label used to identify blank samples
#' @param qc_label Class label for QC sample. If not NULL will use QC samples to calculate the mean intensity
#' @param remove Remove blank samples from peak matrix or not
#' @param fraction_in_blank Number between 0 to 1 to specify fraction in how many blanks peaks should be present.  
#' @export

filter_peaks_by_blank <- function(df, fold_change, classes, blank_label, qc_label=NULL, remove=TRUE, fraction_in_blank=0){
  
  df <- check_peak_matrix_orientation(peak_data = df, classes = classes)
  
  M_blanks <- rbind(df[ ,classes == blank_label], NULL)
  
  if (!is.null(qc_label)){
    M_non_blanks <- df[1:dim(df)[1], classes == qc_label, drop=FALSE]
  } else {
    M_non_blanks <- df[1:dim(df)[1], classes != blank_label, drop=FALSE]
  }

  FUN <- function(x) median(x,na.rm=T)

  median_intensity_blanks <- apply(M_blanks, 1, FUN)
  median_intensity_non_blanks <- apply(M_non_blanks, 1, FUN)

  fold_change_exp <- median_intensity_non_blanks / median_intensity_blanks
  
  blank_fraction <- 1-(apply (is.na(M_blanks), 1, sum)  / ncol(M_blanks))
  
  idxs <- fold_change_exp >= fold_change & blank_fraction >= fraction_in_blank

  idxs[which(is.na(median_intensity_non_blanks))] <- FALSE
  idxs[which(is.na(median_intensity_blanks))] <- TRUE

  if (is.null(qc_label)){
      flags <- cbind(median_non_blanks=median_intensity_non_blanks, median_blanks=median_intensity_blanks, 
                    fold_change=fold_change_exp, blank_flags=as.numeric(idxs), 
                    blank_fraction_flags=as.numeric(blank_fraction))
  } else {
      flags <- cbind(median_QCs=median_intensity_non_blanks, median_blanks=median_intensity_blanks, 
                    fold_change=fold_change_exp, blank_flags=as.numeric(idxs),
                    blank_fraction_flags=as.numeric(blank_fraction))
  }

  df <- df[idxs, ]

  if (remove){
    df <- df[ ,classes != blank_label]
  }
  return(list(df = df, flags = flags))
}

#' Filter features by fraction of missing values
#'
#' @param df Peak intensity matrix
#' @param min_frac Threshold of fraction of detection
#' @param classes Vector of class labels
#' @param method Method to use, default is QC
#' @param qc_label Class label for QC sample
#' @export

filter_peaks_by_fraction = function(df, min_frac, classes=NULL, method="QC", qc_label="QC"){
  
  df <- check_peak_matrix_orientation(peak_data = df, classes = classes)
  FUN <- function(irr) return(length(which(!is.na(irr)))/length(irr))

  if (method == "within" || method == "QC"){
    fracs <- data.frame(matrix(ncol=dim(df)[1], nrow=0))
    fracs_flags <- data.frame(matrix(ncol=dim(df)[1], nrow=0))
    rn_fracs <- c()
    rn_fracs_flags <- c()
    if (method == "within"){
      cls <- unique(classes)
    } else {
      cls <- c(qc_label)
    }
    for (cl in cls){
      idxs <- cl == classes
      subset_cl <- df[ ,idxs]
      frac <- apply(subset_cl, 1, FUN)
      fracs <- rbind(fracs, round(frac,2))
      fracs_flags <- rbind(fracs_flags, as.numeric(frac >= min_frac))
      rn_fracs <- c(rn_fracs, paste("fraction_", cl, sep=""))
      rn_fracs_flags <- c(rn_fracs_flags, paste("fraction_flag_", cl, sep=""))
    }
    colnames(fracs) <- rownames(df)
    colnames(fracs_flags) <- rownames(df)
    rownames(fracs) <- rn_fracs
    rownames(fracs_flags) <- rn_fracs_flags
    flags <- t(rbind(fracs, fracs_flags))
    idxs <- colSums(fracs_flags) > 0
  } else if (method == "across"){
    frac <- apply(df, 1, FUN)
    idxs <- frac >= min_frac
    flags <- cbind(fraction=round(frac, 2), fraction_flags=as.numeric(idxs))
  }
  return(list(df = df[idxs, , drop=FALSE], flags = flags))
}

#' Filter features by RSD% of QC samples
#'
#' @param df Peak intensity matrix
#' @param max_rsd Threshold of QC RSD% value
#' @param classes Vector of class labels
#' @param qc_label Class label for QC sample
#' @export

filter_peaks_by_rsd = function(df, max_rsd, classes, qc_label){
  
  df <- check_peak_matrix_orientation(peak_data = df, classes = classes)
  
  df_qcs <- df[ ,classes == qc_label, drop=FALSE]

  FUN <- function(x) sd(x,na.rm=T)/mean(x,na.rm=T)*100.0
  
  rsd_values <- apply(df_qcs, 1, FUN)

  idxs = rsd_values < max_rsd
  idxs[is.na(idxs)] = FALSE

  flags = cbind(rsd_QC=round(rsd_values,2), rsd_flags=as.numeric(idxs))

  return(list(df = df[idxs, , drop=FALSE], flags = flags))
}

#' Missing values filter
#'
#' @param df Input data frame.
#' @param max_perc_mv Threshold of missing value percentage.
#' @export

filter_samples_by_mv = function(df, max_perc_mv){

  FUN = function(irr) return(length(which(is.na(irr)))/length(irr))
  perc_mv = apply(df,1,FUN)
  idxs =  perc_mv <= max_perc_mv

  flags = cbind(perc_mv=round(perc_mv,2), flags=as.numeric(idxs))

  return(list(df = df[idxs, 1:length(df[1,]), drop=FALSE], flags = flags))
}
