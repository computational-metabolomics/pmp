#' @importFrom stats median sd var 
#' @importFrom S4Vectors DataFrame
#' @import SummarizedExperiment
#' @importFrom matrixStats rowMedians
#' @importFrom matrixStats rowCounts
#' @importFrom matrixStats rowSds

NULL

#' Filter features by blank samples
#'
#' Metabolomics datasets often contain many features of non-biological 
#' origin e.g. those associated with extraction and analysis solvents. 
#' This tool facilitates the removal of such features from the data matrix, 
#' as defined using an appropriate blank sample.
#' 
#' If parameter \code{qc_label} is not \code{NULL}, QC samples which will be 
#' used to calculate the median signal intensity.
#'
#' @param df A matrix-like (e.g. an ordinary matrix, a data frame) or 
#' \link[SummarizedExperiment]{RangedSummarizedExperiment-class} object with 
#' all values of class \code{numeric()} or \code{integer()} of peak 
#' intensities, areas or other quantitative characteristic.
#' @param fold_change \code{numeric(1)}, fold_change minimum fold change 
#' between analytical and blank samples.
#' @param classes \code{character()}, vector of class labels. Must be the same 
#' length as the number of sample in the input peak table. If input is 
#' \code{SummarizedExperiment} object, use 
#' \code{SummarizedExperiment_object$meta_data_column_name}.
#' @param blank_label \code{character(1)}, class label used to identify blank
#' samples.
#' @param qc_label \code{character(1)} or \code{NULL}, class label used to
#' identify QC samples.
#' @param remove_samples \code{logical(1)}, remove blank samples from peak 
#' matrix or not.
#' @param remove_peaks \code{logical(1)}, remove filtered features from peak 
#' matrix or not.
#' @param fraction_in_blank \code{numeric(1)}, value between 0 to 1 to specify 
#' fraction in how many blanks peaks should be present.
#' 
#' @return Object of class \code{SummarizedExperiment}. If input data are a 
#' matrix-like (e.g. an ordinary matrix, a data frame) object, function returns 
#' \code{numeric()} matrix-like object of filtered data set. Function 
#' \code{flags} are added to the object \code{attributes} and is a 
#' \link[S4Vectors]{DataFrame-class} with five columns. The same 
#' \code{DataFrame} object containing flags is added to \code{rowData()} 
#' element of \code{SummarizedExperiment} object as well.\cr
#' \cr
#' Columns in \code{rowData()} or \code{flags} element contain: \cr
#' \code{median_non_blanks} median intensities of features of non-blank
#' samples; \cr
#' \code{median_blanks} median intensities of features of blank samples; \cr
#' \code{fold_change} fold change between analytical and blank samples; \cr
#' \code{blank_flags} \code{integer()}, if 0 feature is flagged to be 
#' removed; \cr
#' \code{blank_fraction_flags} \code{numeric()}, fraction in how many blank
#' samples peaks is present. \cr
#' 
#' @examples
#' df <- MTBLS79[ ,MTBLS79$Batch == 1]
#' df$Class[1:2] <- "Blank"
#' out <- filter_peaks_by_blank(df=df, fold_change=1.2, 
#'    classes=df$Class, blank_label="Blank", qc_label=NULL, 
#'    remove_samples=FALSE, remove_peaks=TRUE, fraction_in_blank=0)
#' 
#' @export

filter_peaks_by_blank <- function(df, fold_change, classes, blank_label, 
        qc_label=NULL, remove_samples=TRUE, remove_peaks=TRUE, 
        fraction_in_blank=0) {
    input_df_class <- class(df)[1]
    df <- check_input_data(df=df, classes=classes)
    M_blanks <- df[, classes == blank_label]
    if (!is.null(qc_label)) {
        M_non_blanks <- df[, classes == qc_label,
            drop=FALSE]
    } else {
        M_non_blanks <- df[, classes != blank_label,
            drop=FALSE]
    }
    median_intensity_blanks <- 
        matrixStats::rowMedians(assay(M_blanks), na.rm=TRUE)
    median_intensity_non_blanks <- 
        matrixStats::rowMedians(assay(M_non_blanks), na.rm=TRUE)
    fold_change_exp <- median_intensity_non_blanks / median_intensity_blanks
    blank_fraction <- 1 - (rowCounts(assay(M_blanks), value=NA) / 
        ncol(assay(M_blanks)))
    idxs <- fold_change_exp >= fold_change & blank_fraction >= fraction_in_blank
    idxs[which(is.na(median_intensity_non_blanks))] <- FALSE
    idxs[which(is.na(median_intensity_blanks))] <- TRUE
    if (is.null(qc_label)) {
        row_data <- DataFrame(median_non_blanks=median_intensity_non_blanks)
    } else {
        row_data <- DataFrame(median_QCs=median_intensity_non_blanks)
    }
    row_data$median_blanks <- median_intensity_blanks
    row_data$fold_change <- fold_change_exp
    row_data$blank_flags <- as.numeric(idxs)
    row_data$blank_fraction_flags <- as.numeric(blank_fraction)
    rowData(df) <- cbind(rowData(df), row_data)
    if (remove_peaks){
        df <- df[idxs, ]
    }
    if (remove_samples) {
        df <- df[, classes != blank_label]
    }
    meta_data <- metadata(df)
    meta_data$processing_history$filter_peaks_by_blank <- return_function_args()
    metadata(df) <- meta_data
    df <- return_original_data_structure(df)
    return(df)
}

#' Filter features by fraction of missing values
#' 
#' Metabolomics datasets often contain 'features' with irreproducible peak 
#' intensity values, or with large numbers of missing values. This tool 
#' facilitates the remove of such features from a data matrix, based upon the 
#' relative proportion (minimum fraction) of samples containing non-missing 
#' values.
#'
#' @inheritParams filter_peaks_by_blank
#' @param min_frac \code{numeric(1)}, value between 0 and 1, a threshold of 
#' fraction of detection.
#' @param method \code{character(1)}, method to use. \code{QC} - within 
#' QC samples, \code{within} - within each sample class or \code{across} - 
#' across all samples.
#'
#' @return Object of class \code{SummarizedExperiment}. If input data are a 
#' matrix-like (e.g. an ordinary matrix, a data frame) object, function returns 
#' \code{numeric()} matrix-like object of filtered data set. Function 
#' \code{flags} are added to the object \code{attributes} and is a 
#' \link[S4Vectors]{DataFrame-class} with five columns. The same 
#' \code{DataFrame} object containing flags is added to \code{rowData()} 
#' element of \code{SummarizedExperiment} object as well.\cr
#' \cr
#' 
#' Columns in \code{rowData()} or \code{flags} element contain \code{fractions} 
#' of missing values per feature within QC samples (mehtod \code{QC}), 
#' across (method \code{across}) or within (mehtod \code{within}) each sample 
#' group.
#' 
#' @examples 
#' df <- MTBLS79[ ,MTBLS79$Batch == 1]
#' out <- filter_peaks_by_fraction(df=df, min_frac=1, 
#'     classes=df$Class, method='QC', qc_label='QC')
#'     
#' out <- filter_peaks_by_fraction(df=df, min_frac=1, 
#'     classes=df$Class, method='across', qc_label='QC')
#' 
#' out <- filter_peaks_by_fraction(df=df, min_frac=1, 
#'     classes=df$Class, method='within', qc_label='QC')
#' 
#' @export

filter_peaks_by_fraction <- function(df, min_frac, classes=NULL,
    method="QC", qc_label="QC", remove_peaks=TRUE) {
    df <- check_input_data(df=df, classes=classes)
    if (method == "within" || method == "QC") {
        fracs <- data.frame(matrix(ncol=dim(df)[1], nrow=0))
        fracs_flags <- data.frame(matrix(ncol=dim(df)[1], nrow=0))
        rn_fracs <- c()
        rn_fracs_flags <- c()
        if (method == "within") {
            cls <- unique(classes)
        } else {
            cls <- c(qc_label)
        }
        for (cl in cls) {
            idxs <- cl == classes
            subset_cl <- df[, idxs, drop=FALSE]
            frac <- rowSums(!is.na(assay(subset_cl)))/ncol(subset_cl)
            fracs <- rbind(fracs, round(frac, 2))
            fracs_flags <- rbind(fracs_flags, as.numeric(frac >= min_frac))
            rn_fracs <- c(rn_fracs, paste("fraction_", cl, sep=""))
            rn_fracs_flags <- c(rn_fracs_flags, paste("fraction_flag_", cl,
                sep=""))
        }
        colnames(fracs) <- rownames(df)
        colnames(fracs_flags) <- rownames(df)
        rownames(fracs) <- rn_fracs
        rownames(fracs_flags) <- rn_fracs_flags
        flags <- t(rbind(fracs, fracs_flags))
        idxs <- colSums(fracs_flags) > 0
    } else if (method == "across") {
        frac <- rowSums(!is.na(assay(df)))/ncol(df)
        idxs <- frac >= min_frac
        flags <- cbind(fraction=round(frac, 2),
            fraction_flags=as.numeric(idxs))
    }
    rowData(df) <- cbind(rowData(df), DataFrame(flags))
    if (remove_peaks){
        df <- df[idxs, ]
    }
    meta_data <- metadata(df)
    function_arguments <- return_function_args()
    filter_name <- paste("filter_peaks_by_fraction", 
        function_arguments$method, sep="_")
    meta_data$processing_history[[filter_name]] <- 
        function_arguments
    metadata(df) <- meta_data
    df <- return_original_data_structure(df)
    return(df)
}

#' Remove features from peak intensity matrix
#'
#' Filter to remove features.
#'
#' @inheritParams filter_peaks_by_blank
#' @param rem_index \code{logical()}, vector containing \code{TRUE} vales for 
#' features to remove. Should be the same length as number of features in input 
#' data.
#' 
#' @return Object of class \code{SummarizedExperiment}. If input data are a 
#' matrix-like (e.g. an ordinary matrix, a data frame) object, function returns 
#' the same R data structure as input with all value of data type 
#' \code{numeric()}.
#' 
#' @examples
#' df <- MTBLS79[ ,MTBLS79$Batch == 1]
#' rem_index <- vector(mode="logical", 
#'     length=nrow(SummarizedExperiment::assay(df)))
#' rem_index[c(1, 20, 456, 789)] <- TRUE
#' out <- remove_peaks(df=df, rem_index=rem_index)
#' 
#' @export

remove_peaks <- function(df, rem_index) {
    df <- check_input_data(df=df)
    if (is.logical(rem_index)) {
        df <- df[!rem_index, ]
        meta_data <- metadata(df)
        meta_data$processing_history$remove_peaks <- 
            return_function_args()
        metadata(df) <- meta_data
        df <- return_original_data_structure(df)
        return(df)
    } else {
        stop("Vector of indexes to remove from peak matrix should be logical
    vector of TRUE/FALSE vaues.")
    }
}

#' Filter features by RSD\% of QC samples
#'
#' Metabolomics datasets often contain 'features' with irreproducible peak 
#' intensity values, or with large numbers of missing values. This tool 
#' facilitates the remove of such features from a data matrix, based upon 
#' relative standard deviation of intensity values for a given feature within 
#' specified QC samples.
#'
#' @inheritParams filter_peaks_by_blank
#' @param max_rsd \code{numeric()}, threshold of QC RSD\% value
#' 
#' @return Object of class \code{SummarizedExperiment}. If input data are a 
#' matrix-like (e.g. an ordinary matrix, a data frame) object, function returns 
#' \code{numeric()} matrix-like object of filtered data set. Function 
#' \code{flags} are added to the object \code{attributes} and is a 
#' \link[S4Vectors]{DataFrame-class} with five columns. The same 
#' \code{DataFrame-class} object containing flags is added to \code{rowData()} 
#' element of \code{SummarizedExperiment} object as well.\cr
#' \cr
#' Columns in \code{rowData()} or \code{flags} element contain: \cr
#' \code{rsd_QC} \code{numeric()}, RSD\% value of QC samples per feature; \cr
#' \code{rsd_flags} \code{integer()},if 0 feature is flagged to be removed. \cr
#' 
#' @examples 
#' 
#' df <- MTBLS79[ ,MTBLS79$Batch == 1]
#' out <- filter_peaks_by_rsd(df=df, max_rsd=20,
#'     classes=df$Class, qc_label='QC')
#' 
#' @export

filter_peaks_by_rsd <- function(df, max_rsd, classes, qc_label,
    remove_peaks=TRUE) {
    df <- check_input_data(df=df, classes=classes)
    df_qcs <- df[, classes == qc_label, drop=FALSE]
    rsd_values <- rowSds(assay(df_qcs), na.rm=TRUE) / 
        rowMeans(assay(df_qcs), na.rm=TRUE) * 100
    idxs <- rsd_values < max_rsd
    idxs[is.na(idxs)] <- FALSE
    row_data <- DataFrame(rsd_QC=round(rsd_values, 2),
        rsd_flags=as.numeric(idxs))
    rowData(df) <- cbind(rowData(df), row_data)
    if (remove_peaks){
        df <- df[idxs, ]
    }
    meta_data <- metadata(df)
    meta_data$processing_history$filter_peaks_by_rsd <- return_function_args()
    metadata(df) <- meta_data
    df <- return_original_data_structure(df)
    return(df)
}

#' Filter samples by missing values
#'
#' Missing values in mass spectrometry metabolomic datasets occur widely and 
#' can originate from a number of sources, including for both technical and 
#' biological reasons. In order for robust conclusions to be drawn from 
#' down-stream statistical testing procedures, the issue of missing values 
#' must first be addressed. This tool facilitates the removal of samples 
#' containing a user-defined maximum percentage of missing values.
#'
#' @inheritParams filter_peaks_by_blank
#' @param max_perc_mv \code{numeric(1)}, Value between 0 and 1 of threshold 
#' of missing value percentage in sample.
#' 
#' @return Object of class \code{SummarizedExperiment}. If input data are a 
#' matrix-like (e.g. an ordinary matrix, a data frame) object, function returns 
#' \code{numeric()} matrix-like object of filtered data set. Function 
#' \code{flags} are added to the object \code{attributes} and is a 
#' \link[S4Vectors]{DataFrame-class} with five columns. The same 
#' \code{DataFrame} object containing flags is added to \code{rowData()} 
#' element of \code{SummarizedExperiment} object as well. If element 
#' \code{colData()} already exists flags are appended to existing values.\cr
#' \cr
#' 
#' Columns in \code{colData()} or \code{flags} element contain: \cr
#' \code{perc_mv} \code{numeric()}, fraction of missing values per sample; \cr
#' \code{flags} \code{integer()},if 0 feature is flagged to be removed. \cr
#' 
#' @examples 
#' 
#' df <- MTBLS79
#' out <- filter_samples_by_mv (df=df, max_perc_mv=0.8)
#' 
#' @export

filter_samples_by_mv <- function(df, max_perc_mv, classes=NULL, 
    remove_samples=TRUE) {
    df <- check_input_data(df=df, classes=classes)
    perc_mv <- colSums(is.na(assay(df)))/nrow(df)
    idxs <- perc_mv <= max_perc_mv
    col_data <- DataFrame(filter_samples_by_mv_perc=round(perc_mv, 2),
        filter_samples_by_mv_flags=as.numeric(idxs))
    colData(df) <- cbind(colData(df), col_data)
    if (remove_samples){
        df <- df[ ,idxs]
    }
    meta_data <- metadata(df)
    meta_data$processing_history$filter_samples_by_mv <- return_function_args()
    metadata(df) <- meta_data
    df <- return_original_data_structure(df)
    if (!is(df, "SummarizedExperiment")){
        attributes(df)$flags <- as.matrix(col_data)
    }
    return(df)
}
