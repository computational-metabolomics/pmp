#' @importFrom matrixStats rowAnyMissings

NULL

#' Normalisation by total sum of the features per sample
#' 
#' For each sample, every feature intensity value is divided by the total sum of
#' all feature intensity values measured in that sample (\code{NA} values
#' ignored by default), before multiplication by 100; the unit is \%.
#' 
#' @inheritParams mv_imputation
#' 
#' @return Object of class \code{SummarizedExperiment}. If input data are a 
#' matrix-like (e.g. an ordinary matrix, a data frame) object, function returns 
#' the same R data structure as input with all value of data type 
#' \code{numeric()}.
#' 
#' @examples 
#' df <- MTBLS79[ ,MTBLS79$Batch == 1]
#' out <- normalise_to_sum (df=df)
#'
#' @export

normalise_to_sum <- function(df, check_df=TRUE) {
    if (check_df == TRUE) {
        df <- check_input_data(df)
    } else {
        # normalise_to_sum doesn't need class labels
        # Create generic class label vector to avoid DF to be transposed
        df <- check_input_data(df=df, classes=rep("S", ncol(df)))
    }
    assay(df) <- (sweep(assay(df), 2, colSums(assay(df), na.rm=TRUE)/100, 
        FUN="/"))
    meta_data <- metadata(df)
    meta_data$processing_history$normalise_to_sum <- 
        return_function_args()
        #list (check_df=check_df)
    metadata(df) <- meta_data
    df <- return_original_data_structure(df)
    df
}

#' Calculate reference mean of samples
#' 
#' @param df_qc \code{numeric()}, peak matrix of QC samples.
#' 
#' @return vector of reference mean values
#' @noRd
calculate_ref_mean <- function(df_qc,ref_method='mean'){
    
    if (ref_method=='mean') {
        ref_mean <- rowMeans(df_qc, na.rm=TRUE)
    } else if (ref_method=='median') {
        ref_mean <- rowMedians(df_qc,na.rm=TRUE)
    } else {
        stop(paste0('Unknown ref_method for calculate_ref_mean: "',ref_method,
            '". Valid options are ',
            '"mean" or "median".'))
    }
    return(ref_mean)
}

#' Probabilistic quotient normalisation (PQN)
#' 
#' For every feature the mean response is calculated across all QC samples. A 
#' reference vector is then generated. The median between the reference vector
#' and every sample is computed obtaining a vector of coefficients related to
#' each sample. Each sample is then divided by the median value of the vector
#' of coefficients; this median value is different for each sample. This 
#' method was adapted by Dieterle et al. (2006) (see references). Its purpose 
#' is to take into account the concentration changes of some metabolite 
#' features that affect limited regions of the data.
#'
#' @references Dieterle F. et al., Anal. Chem., 78(13), 2006. 
#' http://dx.doi.org/10.1021/ac051632c
#'
#' @inheritParams filter_peaks_by_blank 
#' @param ref_mean \code{numeric()} or \code{NULL}, Vector of reference mean
#' values to use instead of calculating from QC sample group. If set to 
#' \code{NULL}, QC sample data will be used.
#' @param qc_frac \code{numeric()} A value between 0 and 1 to indicate the 
#' minimum proportion of QC samples a feature must be present in for it to be included 
#' when computing the reference. Default \code{qc_frac = 0}.
#' @param sample_frac \code{numeric()} A value between 0 and 1 to indicate the 
#' minimum proportion of samples a feature must be present in for it to be considered 
#' when computing the normalisation coefficients. Default \code{sample_frac = 0}.
#' @param ref_method \code{character()} Method used to compute the reference from
#' the QC samples. Default \code{ref_method = 'mean'}. Allowed
#' values are "mean" or "median".
#' @return Object of class \code{SummarizedExperiment}. If input data are a 
#' matrix-like (e.g. an ordinary matrix, a data frame) object, function returns 
#' the same R data structure as input with all value of data type 
#' \code{numeric()}.
#' 
#' @examples 
#' df <- MTBLS79[ , MTBLS79$Batch==1]
#' pqn_normalisation(df=df,
#'     classes=df$Class, qc_label='QC')
#' 
#' @export

pqn_normalisation <- function(df, classes, qc_label, ref_mean=NULL, qc_frac=0,
    sample_frac=0,ref_method='mean') {
    
    # check df for issues
    df <- check_input_data(df=df, classes=classes)
    
    # add some rownames if not present
    rm_rownames=FALSE
    if (is.null(rownames(assay(df)))) {
        rownames(df)=as.character(1:nrow(df))
        rm_rownames=TRUE
    }
    
    # if no reference then create one
    if (is.null(ref_mean)){
        if (qc_label == "all") {
            # use all samples
            ref <- df
            ref_class = classes
        } else {
            # use only samples with the provided QC label
            ref <- df[, classes == qc_label]
            ref_class = classes[classes == qc_label]
        }
        ## apply missing value filter to reference
        # only keep features present in qc_frac*100% samples
        ref = filter_peaks_by_fraction(
            df=ref,
            min_frac = qc_frac,
            classes = ref_class,
            method='across',
            qc_label=qc_label,
            remove_peaks = TRUE
        )
        # check we have some features left
        if (nrow(assay(ref))<1) {
            stop(paste0('QC filtering was too strict and none of the selected ',
            'reference features survived the filtering. Try reducing the value',
            ' of qc_frac.')
            )
        }
        # record the number of reference samples used
        n_ref = nrow(assay(ref))
        
        # average the reference samples
        ref_mean <- calculate_ref_mean(df_qc=assay(ref),ref_method)
    }
    
    # filter the samples before calculating coefficient
    df_filt = filter_peaks_by_fraction(
        df=df,
        min_frac = sample_frac,
        classes = classes,
        method='across',
        qc_label=qc_label,
        remove_peaks = TRUE
    )
    
    # check we have some features left
    if (nrow(assay(df_filt))<1) {
        stop(paste0('Sample filtering was too strict and none of the',
        ' features survived the filtering. Try reducing the ',
        'value of sample_frac.')
        )
    }

    ## only keep features that survive both filters
    U=intersect(rownames(assay(df_filt)),rownames(assay(ref)))
    ref_mean=ref_mean[U]
    df_filt=df_filt[U,]
    
    # check we have some features left
    if (nrow(assay(df_filt))<1) {
        stop(paste0('No features remain after filtering that are common to ',
        'both the reference and the sample matrix. Try again with lower ',
        'values for QC_frac and/or sample_frac.')
        )
    }
    # record number of features
    n_samp=nrow(df_filt)
    
    ## calculate coefficients
    # divide each sample by reference
    coef = apply(assay(df_filt),2,function(x) {
            return(x / ref_mean)
        }
    )
    # count number of coefficients per sample
    coef_count = apply(coef,2,function(x){
        return(sum(!is.na(x)))
    })
    # median coefficient for each sample
    coef = matrixStats::colMedians(coef,na.rm=TRUE)
    # convert to matrix
    coef=matrix(coef,nrow=nrow(df),ncol=length(coef),byrow = TRUE)
    
    # apply normalisation
    assay(df) <- assay(df) / coef 
    col_data <- DataFrame(pqn_coef=coef[1,])
    colData(df) <- cbind(colData(df), col_data)
    meta_data <- metadata(df)
    meta_data$processing_history$pqn_normalisation <- c(
        return_function_args(),
        list('reference_feature_count'=n_ref,
             'sample_feature_count',n_samp,
             'coef_feature_count',coef_count,
             'ref_mean',ref_mean
        )
    )
    metadata(df) <- meta_data
    df <- return_original_data_structure(df)
    if (!is(df, "SummarizedExperiment")){
        attributes(df)$flags <- as.matrix(col_data)
    }
    
    # remove the temporary rownames if we added them earlier
    if (rm_rownames) {
        rownames(df)=NULL
    }
    
    return(df)
}
