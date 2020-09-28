#' @title A function implementing LD score regression.
#'
#' @param dat1: formmated summary statistics for trait 1.
#' @param dat2: formmated summary statistics for trait 2.
#' @param trait1.name: specify the name of trait 1, default `exposure`.
#' @param trait2.name: specify the name of trait 2, default `outcome`.
#' @param h2.fix.intercept: whether to fix LD score regression intercept to 1, default `FALSE`.
#' @param ldscore.dir: specify the path to the LD score files.
#'
#' @return List with the following elements:
#' \describe{
#' \item{Mdat}{Homonised data set }
#' \item{Sigma_err}{the estimated residual covariance matrix of Z scores.}
#' \item{Omega}{the estimated covariance matrix of background effects}
#' }
#'
#' @export
#'
#'
est_paras <- function(dat1,
                      dat2,
                      trait1.name = "exposure",
                      trait2.name = "outcome",
                      LDSC = T,
                      h2.fix.intercept = F,
                      ldscore.dir = NULLL){

  message("Merge dat1 and dat2 by SNP ...")
  dat = merge(dat1, dat2, by="SNP")

  message("Harmonise the direction of SNP effects of trait 1 and trait 2")
  flip.index = which((dat$A1.x == dat$A2.y & dat$A1.y == dat$A2.x) |
                       (dat$A1.x ==comple(dat$A2.y) & dat$A1.y == comple(dat$A2.x)))

  dat[,"A1.y"] = dat[,"A1.x"]
  dat[,"A2.y"] = dat[,"A2.x"]
  dat[flip.index ,"Z.y"] = -dat[flip.index ,"Z.y"]


  message("Read in LD scores ... ")
  if(is.null(ldscore.dir)) stop("Please provide the information on LD scores")

  ld <- suppressMessages(readr::read_delim(paste0(ldscore.dir,"/1.l2.ldscore.gz"), "\t", escape_double = FALSE, trim_ws = TRUE,progress = F))

  for(i in 2:22){
    ld <- rbind(ld,suppressMessages(readr::read_delim(paste0(ldscore.dir, "/", i,".l2.ldscore.gz"), "\t", escape_double = FALSE, trim_ws = TRUE,progress = F)))
  }

  m.chr  <- suppressMessages(readr::read_csv(paste0(ldscore.dir,"/1.l2.M_5_50"),  col_names = FALSE))

  for(i in 2:22){
    m.chr <- rbind(m.chr,suppressMessages(readr::read_csv(paste0(ldscore.dir, "/", i,".l2.M_5_50"),  col_names = FALSE)))
  }

  M = sum(m.chr)  # the number of SNPs include in the LD score estimation

  message("Add LD scores to the harmonised data sets...")
  merged  = merge(dat, ld, by="SNP")

  message("The Harmonised dataset will also be used for  MR analysis \n")
  dat = data.frame(SNP = merged$SNP,
                   A1 = merged$A1.x,
                   A2 = merged$A2.x,
                   b.exp = merged$Z.x/sqrt(merged$N.x),
                   b.out = merged$Z.y/sqrt(merged$N.y),
                   se.exp = 1/sqrt(merged$N.x),
                   se.out = 1/sqrt(merged$N.y),
                   pval.exp = merged$P.x,
                   pval.out = merged$P.y,
                   L2 = merged$L2)
  
  
Sigma_err = NULL
ldsc_res  = NULL
Omega = NULL
  
if(LDSC){
  
  message("Begin estimation of Sigma and Omega using LDSC ...")
  gcres12 = ldsc_GC(merged,
                    trait1.name = "exposure",
                    trait2.name = "outcome",
                    Twostep = T,
                    M=M,
                    h2.fix.intercept = h2.fix.intercept)

  Sigma_err = matrix(as.vector(gcres12$I), nrow=2, ncol=2)
  Omega = matrix(as.vector(gcres12$cov)/M, nrow=2, ncol=2)
  
  }

  return(list(dat=dat,
              ldsc_res=gcres12,
              Sigma_err = Sigma_err,
              Omega = Omega))
}




