#' @title Format GWAS summary data.
#' @description  Reads in GWAS summary data. Infer Zscores from p-values and signed satatistics.
#' This function is adapted from the format_data() function in MRCIEU/TwoSampleMR.
#'
#' @md
#' @param dat Data frame. Must have header with at least SNP A1 A2 signed statistics pvalue and sample size.
#' @param snps.merge Data frame with SNPs to extract. must have headers: SNP A1 and A2. For example, the hapmap3 SNPlist.
#' @param snps.remove a set of SNPs needed to be removed. For example, the SNPs in MHC region.
#' @param snp_col column with SNP rs IDs. The default is `NULL`.
#' @param b_col   Name of column with effect sizes. The default is `NULL`.
#' @param or_col: Name of column with odds ratio. The default is `NULL`.
#' @param se_col Name of column with standard errors. The default is `NULL`.
#' @param freq_col Name of column with effect allele frequency. The default is `NULL`.
#' @param A1_col  Name of column with effect allele. Must contain only the characters "A", "C", "T" or "G". The default is `NULL`.
#' @param A2_col  Name of column with non effect allele. Must contain only the characters "A", "C", "T" or "G". The default is `NULL`.
#' @param p_col  Name of column with p-value. The default is `NULL`.
#' @param ncase_col Name of column with number of cases. The default is `NULL`.
#' @param ncontrol_col Name of column with number of controls. The default is `NULL`.
#' @param n_col Name of column with sample size. The default is `NULL`.
#' @param z_col Name of column with Zscore. The default is `NULL`.
#' @param info_col Name of column with inputation Info. The default is `NULL`.
#' @param log_pval The pval is -log10(p_col). The default is `FALSE`.
#' @param min_freq SNPs with allele frequecy less than min_freq will be removed.The default is `0.05`
#' @param n  Sample size
#' @param chi2_max SNPs with tested chi^2 statistics large than chi2_max will be removed.The default is `80`
#' @param n_qc Whether to remove SNPs according to the sample size of SNPs. The default is `FALSE`.
#'
#' @export
#' @return data frame wih headers: SNP: rsid; A1: effect allele; A2: non effect allel;
#' Z: Z score;  N: sample size;  chi2: chi square statistics;  P: p-value.
#' @importFrom stats pnorm
format_data <- function(dat,
                        snps.merge = w_hm3.snplist,
                        snps.remove = MHC.SNPs,
                        snp_col=NULL,
                        b_col= NULL,
                        or_col =NULL,
                        se_col=NULL,
                        freq_col=NULL,
                        A1_col=NULL,
                        A2_col=NULL,
                        p_col=NULL,
                        ncase_col=NULL,
                        ncontrol_col=NULL,
                        n_col=NULL,
                        n=NULL,
                        z_col=NULL,
                        info_col=NULL,
                        log_pval=FALSE,
                        n_qc=F,
                        chi2_max = 80,
                        min_freq=0.05)
{
  message("Begin formatting .... ")
  message("The raw dataset has ", nrow(dat), " dat lines.")
  cols = c(snp_col, b_col, or_col,  se_col, freq_col, A1_col, A2_col,  p_col,  ncase_col,  ncontrol_col, n_col,z_col,info_col)
  dat = dat[, names(dat) %in% cols]

  if(! snp_col %in% names(dat)){
    stop("SNP column not found")
  }

  if(! A1_col %in% names(dat) | (!A2_col %in% names(dat))){
    stop("A1/A2 columns not found")
  }

  if((!b_col %in% names(dat)) & (!or_col %in% names(dat)) & (!z_col %in% names(dat))){
    stop("signed statistics not found")
  }

  if((!n_col %in% names(dat)) & (!(ncase_col %in% names(dat) & ncontrol_col %in% names(dat))) & is.null(n)){
    stop("Information for sample size not found")
  }

  names(dat)[which(names(dat) == snp_col)[1]] <- "SNP"
  dat$SNP <- tolower(dat$SNP)
  dat$SNP <- gsub("[[:space:]]", "", dat$SNP)
  dat <- subset(dat, !is.na(SNP))

  if(!is.null(snps.remove)){
    message("Remove SNPs in MHC region ...")
    dat <- subset(dat, !SNP %in% snps.remove)
  }

  if(!is.null(snps.merge)){
    snps.merge = snps.merge[, c("SNP","A1","A2")]
    colnames(snps.merge) = c("SNP","ref.A1", "ref.A2")
    message("Merge SNPs with the hapmap3 snplist ...")
    dat <- merge(dat, snps.merge, by="SNP")
  }


  if(info_col %in% names(dat)){
    names(dat)[which(names(dat) == info_col)[1]] <- "info"
    message("Remove SNPs with imputation info less than 0.9 ...")
    dat <- subset(dat, info > 0.9)
  }

  # Check effect_allele (A1)
  if(A1_col %in% names(dat)){

    names(dat)[which(names(dat) == A1_col)[1]] <- "A1"

    if(is.logical(dat$A1)){
      dat$A1 <- substr(as.character(dat$A1), 1, 1)
    }

    if(!is.character(dat$A1)){
      message("effect_allele column is not character data. Coercing...")
      dat$A1 <- as.character(dat$A1)
    }

    dat$A1 <- toupper(dat$A1)

    index = !grepl("^[ACTG]+$", dat$A1)
    if(any(index)){
      message("effect_allele column has some values that are not A/C/T/G. Remove these SNPs...")
      dat = dat[-index,]
    }

  }

  # Check effect_allele (A2)
  if(A2_col %in% names(dat)){

    names(dat)[which(names(dat) == A2_col)[1]] <- "A2"

    if(is.logical(dat$A2))
    {
      dat$A2 <- substr(as.character(dat$A2), 1, 1)
    }
    if(!is.character(dat$A2))
    {
      message("other_allele column is not character data. Coercing...")
      dat$A2 <- as.character(dat$A2)
    }

    dat$A2 <- toupper(dat$A2)

    index = !grepl("^[ACTG]+$", dat$A2)
    if(any(index)){
      message("other allele column has some values that are not A/C/T/G. Remove these SNPs...")
      dat = dat[-index,]
    }
  }

  index = which((dat$A1=="A" & dat$A2=="T") |
                  (dat$A1=="T" & dat$A2=="A") |
                  (dat$A1=="C" & dat$A2=="G") |
                  (dat$A1=="G" & dat$A2=="C") |
                  (dat$A1=="A" & dat$A2=="A") |
                  (dat$A1=="T" & dat$A2=="T") |
                  (dat$A1=="C" & dat$A2=="C") |
                  (dat$A1=="G" & dat$A2=="G"))
  if(any(index)){
    message("Remove ambiguous SNPs ...")
    dat = dat[-index,]
    rm(index)
  }

  comple <- function(allele){

    ifelse(allele == "A","T", ifelse(allele == "T","A", ifelse(allele == "G","C", ifelse(allele == "C","G", allele)) ))

  }

  index = which(!((dat$A1 == dat$ref.A2 & dat$A2 == dat$ref.A1) |
                    (dat$A1 == dat$ref.A1 & dat$A2 == dat$ref.A2) |
                    (dat$A1 == comple(dat$ref.A2) & dat$A2 == comple(dat$ref.A1))|
                    (dat$A1 == comple(dat$ref.A1) & dat$A2 == comple(dat$ref.A2))))

  if(any(index)){
    message("Remove SNPs with alleles not matched with the hapmap3 snplist")
    dat = dat[-index,]
    rm(index)
  }

  # message("Remove Duplicated SNPs if have. Just keeping the first instance")
  # dat = dat[unique(dat$SNP),]  # can be faster

  # Check effect size estimate (b)
  # set b as log(or) of b is not available
  if(! b_col %in% names(dat) & or_col %in% names(dat)){
     names(dat)[which(names(dat) == or_col)[1]] <- "or"
    message("infer b column from log(or)...")
    dat$b = log(dat$or)
    
      if(se_col %in% names(dat)){
    names(dat)[which(names(dat) == se_col)[1]] <- "se"
    if(!is.numeric(dat$se))
    {
      message("se column is not numeric. Coercing...")
      dat$se <- as.numeric(dat$se)
    }
    dat = dat[is.finite(dat$se) & dat$se > 0,]
    dat$se <- log(dat$se)
  }
    
  }

  if(b_col %in% names(dat)){
    names(dat)[which(names(dat) == b_col)[1]] <- "b"
    if(!is.numeric(dat$b))
    {
      message("b column is not numeric. Coercing...")
      dat$b <- as.numeric(dat$b)
    }
    dat = dat[is.finite(dat$b),]
      # Check se
  if(se_col %in% names(dat)){
    names(dat)[which(names(dat) == se_col)[1]] <- "se"
    if(!is.numeric(dat$se))
    {
      message("se column is not numeric. Coercing...")
      dat$se <- as.numeric(dat$se)
    }
    dat = dat[is.finite(dat$se) & dat$se > 0,]
  }

  }


  if(z_col %in% names(dat)){
    names(dat)[which(names(dat) == z_col)[1]] <- "Z"
  }

  # Check freq
  if(freq_col %in% names(dat)){
    names(dat)[which(names(dat) == freq_col)[1]] <- "freq"
    if(!is.numeric(dat$freq))
    {
      message("freq column is not numeric. Coercing...")
      dat$freq <- as.numeric(dat$freq)
    }
    # remove SNP with allele freqcy less than min_freq
    dat = dat[dat$freq > min_freq & dat$freq < 1-min_freq, ]
  }


  if(ncase_col %in% names(dat)){
    names(dat)[which(names(dat) == ncase_col)[1]] <- "ncase"
    if(!is.numeric(dat$ncase))
    {
      message(ncase_col, " column is not numeric")
      dat$ncase <- as.numeric(dat$ncase)
    }
  }

  if(ncontrol_col %in% names(dat)){
    names(dat)[which(names(dat) == ncontrol_col)[1]] <- "ncontrol"
    if(!is.numeric(dat$ncontrol))
    {
      message(ncontrol_col, " column is not numeric")
      dat$ncontrol <- as.numeric(dat$ncontrol)
    }
  }

  if(n_col %in% names(dat)){
    names(dat)[which(names(dat) == n_col)[1]] <- "N"

    if(!is.numeric(dat$N)){
      message(samplesize_col, " column is not numeric")
      dat$N <- as.numeric(dat$N)
    }

    if("ncontrol" %in% names(dat) & "ncase" %in% names(dat)){
      index <- is.na(dat$N) & (!is.na(dat$ncase)) & (!is.na(dat$ncontrol))
      if(any(index)){
        message("Generating sample size from ncase and ncontrol")
        dat$N[index] <- dat$ncase[index] + dat$ncontrol[index]
      }
    }
  }else if("ncontrol" %in% names(dat) & "ncase" %in% names(dat)){
    message("Generating sample size from ncase and ncontrol")
    dat$N <- dat$ncase + dat$ncontrol
  }

  if(!"N" %in% names(dat) & (!is.null(N))){
    message("Generating sample size from specified sample size")
    dat$N = N
  }


  # Check pval
  if(p_col %in% names(dat) & log_pval){
    names(dat)[which(names(dat) == p_col)[1]] <- "P"
    dat$P <- 10^-dat$P
  }

  if(p_col %in% names(dat)){
    names(dat)[which(names(dat) == p_col)[1]] <- "P"
    if(!is.numeric(dat$P)){
      message("pval column is not numeric. Coercing...")
      dat$P <- as.numeric(dat$P)
    }
    message("Remove SNPs with P value < 0 or P value > 1")
    dat = subset(dat, P>=0 & P <=1)
  }

  # change colnames to upcase
  dat <- setNames(dat, toupper(names(dat)))

  if("z" %in% names(dat)){
    dat$chi2 = dat$Z^2
  }

    # calculate Z from p value
  if("p" %in% names(dat)){
    if("b" %in% names(dat) & ! "z" %in% names(dat)){
      dat$chi2 = qchisq(dat$P,1,lower.tail = F)
      message("Infer z score from P value and b ...")
      dat$Z = sign(dat$b)* sqrt(dat$chi2)
    }
  }

   # calculate z if not contain z, but with b se or p and sign
  if((! "Z" %in% names(dat))  & ("b" %in% names(dat) & "se" %in% names(dat))){
    message("Infer Z score from b/se ...")
    dat$Z = dat$b/dat$se
    dat$chi2 = dat$Z^2
  }


  # calculate chi2 if no chi2
  if(!"chi2" %in% names(dat)) dat$chi2 = dat$Z^2

  # calculate P if not contain P
  if(!"P" %in% names(dat)) dat$P = pchisq(dat$chi2, 1, lower.tail = F)

  if(! "Z" %in% names(dat)){
    stop("Error: No information for Zscore ")
  }

  if(n_qc==T){
    n_min = quantile(dat$N, 0.9)/1.5
    message("Remove SNPs with sample size <= n_min... ")
    dat = subset(dat, N > n_min)
  }


  message("Remove SNPs with chi2 > chi2_max ... ")
  dat = subset(dat, chi2 < chi2_max)

  message("The formatted data has ", nrow(dat), " dat lines. \n")
  
  dat = dat[, c("SNP","A1","A2","Z","N","chi2","P")]
  
  return(dat)

}


