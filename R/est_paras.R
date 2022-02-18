#' @title A function harmonizing datasets and estimate background parameters by LD score regression.
#'
#' @param dat1: formatted GWAS summary-level data for trait 1.
#' @param dat2: formatted GWAS summary-level data for trait 2.
#' @param trait1.name: specify the name of trait 1, default `exposure`.
#' @param trait2.name: specify the name of trait 2, default `outcome`.
#' @param LDSC: logical, whether to run LD score regression, default `TRUE`. If `FALSE`, the function will not give the parameter estimates but will do harmonizing.
#' @param h2.fix.intercept: logical, whether to fix LD score regression intercept to 1, default `FALSE`.
#' @param ldscore.dir: specify the path to the LD score files.
#'
#' @return List with the following elements:
#' \describe{
#' \item{Mdat}{Homonised data set }
#' \item{C}{the estimated C matrix capturing the effects of sample structure }
#' \item{Omega}{the estimated variance-covariance matrix for polygenic effects}
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
                      ldscore.dir = NULL，
                      ld=NULL,
                      M=NULL){
  
  dat1 %<>% dplyr::mutate_if(is.integer, as.numeric)
  dat2 %<>% dplyr::mutate_if(is.integer, as.numeric)
  dat1 %<>% dplyr::mutate_if(is.factor, as.character)
  dat2 %<>% dplyr::mutate_if(is.factor, as.character)
      
  message("Merge dat1 and dat2 by SNP ...")
  dat = merge(dat1, dat2, by="SNP")

  message("Harmonize the direction of SNP effects of exposure and outcome")
  flip.index = which((dat$A1.x == dat$A2.y & dat$A1.y == dat$A2.x) |
                       (dat$A1.x ==comple(dat$A2.y) & dat$A1.y == comple(dat$A2.x)))

  dat[,"A1.y"] = dat[,"A1.x"]
  dat[,"A2.y"] = dat[,"A2.x"]
  dat[flip.index ,"Z.y"] = -dat[flip.index ,"Z.y"]


  message("Read in LD scores ... ")
  if(is.null(ldscore.dir)&is.null(ld)) stop("Please provide the information on LD scores")

  if(is.null(ld) & !is.null(ldscore.dir)){
    
  ld <- suppressMessages(readr::read_delim(paste0(ldscore.dir,"/1.l2.ldscore.gz"), "\t", escape_double = FALSE, trim_ws = TRUE,progress = F))

  for(i in 2:22){
    ld <- rbind(ld,suppressMessages(readr::read_delim(paste0(ldscore.dir, "/", i,".l2.ldscore.gz"), "\t", escape_double = FALSE, trim_ws = TRUE,progress = F)))
  }
 }
  
  if(is.null(M) & !is.null(ldscore.dir)){
  m.chr  <- suppressMessages(readr::read_csv(paste0(ldscore.dir,"/1.l2.M_5_50"),  col_names = FALSE))

  for(i in 2:22){
    m.chr <- rbind(m.chr,suppressMessages(readr::read_csv(paste0(ldscore.dir, "/", i,".l2.M_5_50"),  col_names = FALSE)))
  }

  M = sum(m.chr)  # the number of SNPs include in the LD score estimation
 }
  message("Add LD scores to the harmonized data set...")
  merged  = merge(dat, ld, by="SNP")

  message("The Harmonized data set will also be used for MR analysis \n")
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


C = NULL
gcres12  = NULL
Omega = NULL

if(LDSC){

  message("Begin estimation of C and Omega using LDSC ...")
  gcres12 = ldsc_GC(merged,
                    trait1.name = "exposure",
                    trait2.name = "outcome",
                    Twostep = T,
                    M=M,
                    h2.fix.intercept = h2.fix.intercept)

  C = matrix(as.vector(gcres12$I), nrow=2, ncol=2)
  Omega = matrix(as.vector(gcres12$cov)/M, nrow=2, ncol=2)

  }

  return(list(dat=dat,
              ldsc_res=gcres12,
              C = C,
              Omega = Omega))
}




