#'@title  A function for causal effect estimation and inference in sensitivity analysis.
#' @param MRdat  data.frame at least contains the following columns: b.exp b.out se.exp se.out L2 Threshold. L2:LD score, Threshold: modified IV selection threshold for correction of selection bias
#' @param exposure  exposure name
#' @param outcome   outcome name
#' @param pi0 initial value for pi0, default `NULL` will use the default initialize procedure.
#' @param r  The value of r (correlation between IV strength and direct effect in the foreground model), default `0`, r is fixed.
#' @param Sigma initial value for Sigma (the variance-covariance matrix for forground effects), default `NULL`will use the default initialize procedure.
#' @param C  the estimated C matrix capturing the effects of sample structure. default `diag(2)`.
#' @param Omega  the estimated variance-covariance matrix of polygenic effects. default `matrix(0,2,2)`.
#' @param tol     tolerence, default '1e-12'.
#' @param Cor.SelectionBias   logical, whether to correct selection bias or not. If FALSE, the model won't correct for selection bias.

#' @return a list with the following elements:
#' \describe{
#' \item{exposure: }{exposure of interest}
#' \item{outcome: }{outcome of interest}
#' \item{beta: }{causal effect estimate}
#' \item{beta.se: }{standard error}
#' \item{pval: }{p-value}
#' }
#'
#' @examples
#' library(MRAPSS)
#' exposure = "BMI"
#' outcome = "T2D"
#' Threshold = 5e-05  # IV selection Threshold
#' data(C)
#' data(Omega)
#' data(MRdat)
#' MRres = MR_fixr(MRdat,
#'                 exposure = "BMI",
#'                 outcome = "T2D",
#'                 r=0,
#'                 C = C,
#'                 Omega =  Omega ,
#'                 Cor.SelectionBias = T)
#' @export

MR_fixr <- function(MRdat = NULL,
                        exposure = "exposure",
                        outcome = "outcome",
                        pi0 = NULL,
                        Sigma=NULL,
                        r=0,
                        C = matrix(c(1,0,0,1), 2, 2),
                        Omega = matrix(0, 2, 2),
                        Cor.SelectionBias = T,
                        tol = 1e-08){

  if(is.null(MRdat)){
    cat("No data for MR testing")
    return(NULL)
  }

  if(nrow(MRdat) < 4) stop(" Not enough IVs.")


  if(!Cor.SelectionBias){

    Threshold = 1
    message("Threshold = 1, the model will not correct for selection bias")

  }else{

    Threshold = unique(MRdat$Threshold)

    if(is.null(Threshold)) Threshold = max(MRdat$pval.exp)

  }


  m = nrow(MRdat)

  ## stage 1
  fit_s1 = MR_EM_fixr_func(MRdat,
                               fix.beta = T,
                               beta = 0,
                               r=r,
                               pi0 = pi0,
                               Sigma=Sigma,
                               C = C,
                               Omega = Omega,
                               tol = tol,
                               Threshold = Threshold)

  # stage 2
  fit_s2 = MR_EM_fixr_func(MRdat,
                           fix.beta = F,
                           beta = 0,
                               pi0 = fit_s1$pi0,
                               r=r,
                               Sigma = fit_s1$Sigma,
                               C = C,
                               Omega = Omega,
                               tol =  tol,
                               Threshold = Threshold)
  # Inference
  LR = 2*(fit_s2$likeli-fit_s1$likeli)
  pvalue = pchisq(LR,1,lower.tail = F)
  pvalue = formatC(pvalue, format = "e", digits = 4)
  beta.se = suppressWarnings(abs(fit_s2$beta/sqrt(LR)))


  #cat("***********************************************************\n")
  #cat("MR test results of ", exposure , " on ", outcome, ": \n")
  #cat("r=", r, "beta = ", round(fit_s2$beta,4), "beta.se = ", round(beta.se, 4), "pvalue = ", pvalue, "#SNPs= ", nrow(MRdat), "\n")
  #cat("# valid IVs with foreground signals: ", fit_s2$pi0 * nrow(MRdat), "\n")
  #cat("***********************************************************\n")
  
  return( list(exposure = exposure,
               outcome = outcome,
               beta = fit_s2$beta,
               beta.se = beta.se,
               pvalue = pvalue))
}


