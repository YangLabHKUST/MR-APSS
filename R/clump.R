#' @title Perform LD clumping
#'
#' @description  Peform LD clumping, to prune SNPs in LD within a window. Keep the most significant ones.
#'
#' @param dat  a data frame must have columns with information about SNPs and p values
#' @param SNP_col column with SNP rsid. The default is `"SNP"`
#' @param pval_col column with p value. The default is `"pval"`
#' @param clump_kb clumping window in kb. Default is 1000.
#' @param clump_r2 clumping r2 threshold. Default is 0.001.
#' @param clump_p  clumping significance level for index variants. Default = 5e-05
#' @param bfile     bfile as LD reference panel. If this is provided, then will use local PLINK. Default = NULL.
#' @param plink_bin path to local plink binary. Default = NULL.
#'
#' @return data frame of clumped SNPs
#' @export

clump <- function(dat,
                  IV.Threshold = 5e-05,
                  SNP_col = "SNP",
                  pval_col = "pval.exp",
                  clump_kb = 1000,
                  clump_r2 = 0.001,
                  clump_p = 0.999,
                  pop = "EUR",
                  bfile = NULL,
                  plink_bin = NULL){
  
  dat = dat[which(dat[, pval_col] <= IV.Threshold), ]
    
  df <- data.frame(rsid = dat[, SNP_col], pval = dat[,pval_col])
  colnames(df) = c("rsid", "pval")

  out <- ieugwasr::ld_clump(df, clump_kb=clump_kb, clump_r2=clump_r2, clump_p=clump_p, bfile=bfile, plink_bin = plink_bin, pop=pop)

  MRdat <- dat[which(df$rsid %in% out$rsid),]
  
  if( IV.Threshold <= 5e-07){
    
    MRdat$Threshold =  IV.Threshold
    
  }else{
    
    ratio = ifelse(median(MRdat[, pval_col])/median(dat[, pval_col]) > 1, 1, median(MRdat[, pval_col])/median(dat[, pval_col]))
    
    MRdat$Threshold = ratio * IV.Threshold
    
  }

  return(MRdat)
}
