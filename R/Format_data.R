#' @title Format GWAS summary data.
#' @description  Reads in GWAS summary data. Infer Zscores from p-values and signed satatistics.
#' This function is adapted from the format_data() function in MRCIEU/TwoSampleMR.
#'
#' @md
#' @param dat Data frame. Must have header with at least SNP A1 A2 signed statistics pvalue and sample size.
#' @param snps.merge Data frame with SNPs to extract. must have headers: SNP A1 and A2. For example, the hapmap3 SNPlist.
#' @param snps.remove a set of SNPs needed to be removed. For example, the SNPs in MHC region.
#' @param snp_col column with SNP rs IDs. The default is `SNP`.
#' @param b_col   Name of column with effect sizes. The default is `b`.
#' @param or_col: Name of column with odds ratio. The default is `or`.
#' @param se_col Name of column with standard errors. The default is `se`.
#' @param freq_col Name of column with effect allele frequency. The default is `frew`.
#' @param A1_col  Name of column with effect allele. Must contain only the characters "A", "C", "T" or "G". The default is `A1`.
#' @param A2_col  Name of column with non effect allele. Must contain only the characters "A", "C", "T" or "G". The default is `A2`.
#' @param p_col  Name of column with p-value. The default is `p`.
#' @param ncase_col Name of column with number of cases. The default is `ncase`.
#' @param ncontrol_col Name of column with number of controls. The default is `ncontrol`.
#' @param n_col Name of column with sample size. The default is `n`.
#' @param z_col Name of column with Zscore. The default is `z`.
#' @param info_col Name of column with inputation Info. The default is `info`.
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
                        snp_col="SNP",
                        b_col= "b",
                        or_col = "or",
                        se_col= "se",
                        freq_col= "freq",
                        A1_col="A1",
                        A2_col="A2",
                        p_col="p",
                        ncase_col="ncase",
                        ncontrol_col="ncontrol",
                        n_col="n",
                        n=NULL,
                        z_col="z",
                        info_col="INFO",
                        log_pval=FALSE,
                        chi2_max = NULL,
                        min_freq=0.05){
  
  message("Begin formatting .... ")
  message("The raw dataset has ", nrow(dat), " dat lines")
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

  
  ## Remove NA
  names(dat)[which(names(dat) == snp_col)[1]] <- "SNP"
  dat$SNP <- tolower(dat$SNP)
  dat$SNP <- gsub("[[:space:]]", "", dat$SNP)
  dat <- subset(dat, !is.na(SNP))

  
  ## check info
  if(info_col %in% names(dat)){
    names(dat)[which(names(dat) == info_col)[1]] <- "info"
    dat <- subset(dat, info > 0.9)
    message("Remove SNPs with imputation info less than 0.9 ...", ", remaining ", nrow(dat), " SNPs.")
  }

  
  ## Check effect_allele (A1)
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
       dat = dat[-index,]
       message("effect_allele column has some values that are not A/C/T/G. Remove these SNPs...", ", remaining ", nrow(dat), " SNPs.")
      
    }

  }

  
  ## Check other_allele (A2)
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
      dat = dat[-index,]
      message("other allele column has some values that are not A/C/T/G. Remove these SNPs...", ", remaining ", nrow(dat), " SNPs.")
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
    dat = dat[-index,]
    message("Remove ambiguous SNPs ...", ", remaining ", nrow(dat), " SNPs.")
    rm(index)
  }

  
  ## remove MHC SNPs
   if(!is.null(snps.remove)){
      dat <- subset(dat, !SNP %in% snps.remove)
      message("Remove SNPs in MHC region ...", ", remaining ", nrow(dat), " SNPs.")
  }
  
  
  ## Duplicated SNPs
  dup_SNPs = dat$SNP[which(duplicated(dat$SNP))]
  dat = subset(dat, !SNP %in% dup_SNPs)
  message("Remove  duplicated SNPs ...", ", remaining ", nrow(dat), " SNPs.")

  
  ## merge with hapmap3 SNPlists
  if(!is.null(snps.merge)){
    snps.merge = snps.merge[, c("SNP","A1","A2")]
    colnames(snps.merge) = c("SNP","ref.A1", "ref.A2")
    dat <- merge(dat, snps.merge, by="SNP")
    message("Merge SNPs with the hapmap3 snplist ...", ", remaining ", nrow(dat), " SNPs.")
      comple <- function(allele){
                     ifelse(allele == "A","T", ifelse(allele == "T","A", ifelse(allele == "G","C", ifelse(allele == "C","G", allele)) ))
                 }

    index = which(!((dat$A1 == dat$ref.A2 & dat$A2 == dat$ref.A1) |
                    (dat$A1 == dat$ref.A1 & dat$A2 == dat$ref.A2) |
                    (dat$A1 == comple(dat$ref.A2) & dat$A2 == comple(dat$ref.A1))|
                    (dat$A1 == comple(dat$ref.A1) & dat$A2 == comple(dat$ref.A2))))

    if(any(index)){
      dat = dat[-index,]
      message("Remove SNPs with alleles not matched with the hapmap3 snplist", ", remaining ", nrow(dat), " SNPs.")
       rm(index)
     }
  }
  
  
  ## Check effect size estimate (b) and se
    if(b_col %in% names(dat)){
       names(dat)[which(names(dat) == b_col)[1]] <- "b"
       if(!is.numeric(dat$b)){
          message("b column is not numeric. Coercing...")
          dat$b <- as.numeric(as.character(dat$b))
          }
        dat = dat[is.finite(dat$b),]
        # Check se
        if(se_col %in% names(dat)){
           names(dat)[which(names(dat) == se_col)[1]] <- "se"
           if(!is.numeric(dat$se)){
              message("se column is not numeric. Coercing...")
              dat$se <- as.numeric(as.character(dat$se))
           }
         dat = dat[is.finite(dat$se) & dat$se > 0,]
         }
      }
  
  
  ## set b as log(or) of b is not available
  if(! b_col %in% names(dat) & or_col %in% names(dat)){
     names(dat)[which(names(dat) == or_col)[1]] <- "or"
     message("infer b column from log(or)...")
     dat$b = log(dat$or)
     if(se_col %in% names(dat)){
        names(dat)[which(names(dat) == se_col)[1]] <- "se"
        if(!is.numeric(dat$se)){
          message("se column is not numeric. Coercing...")
          dat$se = as.numeric(as.character(dat$se))
        }
        dat = dat[is.finite(dat$se) & dat$se > 0,]
        dat$se = log(dat$se)
      }
    }

  
  ## Check freq
  if(freq_col %in% names(dat)){
    names(dat)[which(names(dat) == freq_col)[1]] <- "freq"
    if(!is.numeric(dat$freq))
    {
      message("freq column is not numeric. Coercing...")
      dat$freq <- as.numeric(as.character(dat$freq))
    }
    # remove SNP with allele freqcy less than min_freq
    dat = dat[dat$freq > min_freq & dat$freq < 1-min_freq, ]
  }

  
  ## Check n
  if(ncase_col %in% names(dat)){
    names(dat)[which(names(dat) == ncase_col)[1]] <- "ncase"
    if(!is.numeric(dat$ncase))
    {
      message(ncase_col, " column is not numeric")
      dat$ncase <- as.numeric(as.character(dat$ncase))
    }
  }

  if(ncontrol_col %in% names(dat)){
    names(dat)[which(names(dat) == ncontrol_col)[1]] <- "ncontrol"
    if(!is.numeric(dat$ncontrol))
    {
      message(ncontrol_col, " column is not numeric")
      dat$ncontrol <- as.numeric(as.character(dat$ncontrol))
    }
  }

  if(n_col %in% names(dat)){
    names(dat)[which(names(dat) == n_col)[1]] <- "n"

    if(!is.numeric(dat$n)){
      message(samplesize_col, " column is not numeric")
      dat$n <- as.numeric(as.character(dat$n))
    }

    if("ncontrol" %in% names(dat) & "ncase" %in% names(dat)){
      index <- is.na(dat$n) & (!is.na(dat$ncase)) & (!is.na(dat$ncontrol))
      if(any(index)){
        message("Generating sample size from ncase and ncontrol")
        dat$n[index] <- dat$ncase[index] + dat$ncontrol[index]
      }
    }
  }else if("ncontrol" %in% names(dat) & "ncase" %in% names(dat)){
    message("Generating sample size from ncase and ncontrol")
    dat$n <- dat$ncase + dat$ncontrol
  }

  if(!"n" %in% names(dat) & (!is.null(n))){
    message("Generating sample size from specified sample size")
    dat$n = n
  }


  ## Check pval
  if(p_col %in% names(dat) & log_pval){
    names(dat)[which(names(dat) == p_col)[1]] <- "p"
    dat$p <- 10^-dat$p
  }

  if(p_col %in% names(dat)){
    names(dat)[which(names(dat) == p_col)[1]] <- "p"
    if(!is.numeric(dat$p)){
      message("pval column is not numeric. Coercing...")
      dat$p <- as.numeric(dat$p)
    }
    dat = subset(dat, p>=0 & p <=1)
    message("Remove SNPs with p value < 0 or p value > 1", ", remaining ", nrow(dat), " SNPs.")
  }


    
  ## Check z 
  if(z_col %in% names(dat)){
    names(dat)[which(names(dat) == z_col)[1]] <- "z"
  }

  if("z" %in% names(dat)){
    dat$chi2 = dat$z^2
  }
    
  # calculate z from p value
  if("p" %in% names(dat)){
    if("b" %in% names(dat) & ! "z" %in% names(dat)){
      dat$chi2 = qchisq(dat$p,1,lower.tail = F)
      message("Infer z score from p value and b ...")
      dat$z = sign(dat$b)* sqrt(dat$chi2)
    }
  }

  # calculate z if not contain z, but with b se 
  if((! "z" %in% names(dat))  & ("b" %in% names(dat) & "se" %in% names(dat))){
    message("Infer z score from b/se ...")
    dat$z = dat$b/dat$se
    dat$chi2 = dat$z^2
  }

  # calculate chi2 if no chi2
  if(!"chi2" %in% names(dat)) dat$chi2 = dat$z^2

  # calculate p if not contain p
  if(!"p" %in% names(dat)) dat$p = pchisq(dat$chi2, 1, lower.tail = F)

  if(! "z" %in% names(dat)){
    stop("Error: No information for z score ")
  }

  # check missing
  dat = dat[, c("SNP","A1","A2","z","n","chi2","p")] 
  dat = na.omit(dat)
  message("Remove missing values", ", remaining ", nrow(dat), " SNPs.")
  
  n_min = mean(dat$n) - 5* sd(dat$n)
  n_max = mean(dat$n) + 5* sd(dat$n)
  dat = subset(dat, n >= n_min  & n <= n_max)
  message("Remove SNPs with sample size 5 standard deviations away from the mean", ", remaining ", nrow(dat), " SNPs.")
  
  if(is.null(chi2_max)) chi2_max = max(c(80, median(dat$n)/1000))
  dat = subset(dat, chi2 < chi2_max)
  message("Remove SNPs with chi2 > chi2_max ... ", ", remaining ", nrow(dat), " SNPs.")
  
  message("The formatted data has ", nrow(dat), " dat lines. \n")
  

  colnames(dat) = c("SNP","A1","A2","Z","N","chi2","P")
                
  dat %>% dplyr::mutate_if(is.integer, as.numeric)
                
  return(dat)

}


