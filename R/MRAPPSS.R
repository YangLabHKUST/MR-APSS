#'@title  A function for implementing MR-APPSS.
#'@description MR-APPSS: a unified approach to Mendelian Randomization accounting for polygenicity, pleiotropy and sample structure using genome-wide summary statistics.
#'MA-APPSS uses a variantional EM algorithm for estimation of parameters.
#' MR-APPSS uses likelihood ratio test for inference.
#'
#' @param MRdat  data frame at least contain the following varaibles: b.exp b.out se.exp se.out L2 Threshold. L2:LD score, Threshold: modified IV selection threshold for correction of selection bias
#' @param exposure exposure name
#' @param outcome   outcome name
#' @param pi0 initial value for pi0, default `NULL` will use the default initialize procedure.
#' @param sigma.sq initial value for sigma.sq , default `NULL`will use the default initialize procedure.
#' @param tau.sq initial value for tau.sq , default `NULL` will use the default initialize procedure.
#' @param Sigma_err the error term correlation matrix. default `diag(2)`.
#' @param Omega  the background varaince component. default `matrix(0,2,2)`.
#' @param tol     tolerence, default '1e-08'
#' @param Cor.SelectionBias   Whether use the selection Threshold for correction of selection bias. If FALSE, the model won't correct for selection bias.
#' @param ELBO     Whether check the evidence lower bound or not, if `FALSE`, check the maximum likelihood instead. default `FALSE`.
#'
#' @return a list with the following elements:
#' \describe{
#' \item{MRdat: }{Input data frame}
#' \item{exposure: }{exposure of interest}
#' \item{outcome: }{outcome of interest}
#' \item{beta: }{causal effect estimate}
#' \item{beta.se: }{standard error}
#' \item{pval: }{p-value}
#' \item{sigma.sq: }{variance of forground exposure effect}
#' \item{tau.sq: }{variance of forground outcome effect}
#' \item{pi0: }{The probability of a SNP with forground signal after selection}
#' \item{post: }{Posterior estimates of latent varaibles}
#' \item{method: }{"MR-APPSS"}
#' }
#'
#' @examples
#' library(MRAPPSS)
#' exposure = "BMI"
#' outcome = "T2D"
#' Threshold = 5e-05  # IV selection Threshold
#' data(Sigma_err)
#' data(Omega)
#' data(MRdat)
#' MRres = MRAPSS(MRdat,
#'                exposure = "BMI",
#'                outcome = "T2D",
#'                Sigma_err = Sigma_err,
#'                Omega =  Omega ,
#'                Cor.SelectionBias = T)
#' MRplot(MRres, exposure = "BMI", outcome = "T2D")
#' @export

MRAPPSS <- function(MRdat = NULL,
                   exposure = "exposure",
                   outcome = "outcome",
                   pi0 = NULL,
                   sigma.sq = NULL,
                   tau.sq = NULL,
                   Sigma_err = matrix(c(1,0,0,1), 2, 2),
                   Omega = matrix(0, 2, 2),
                   Cor.SelectionBias = T,
                   tol = 1e-08,
                   ELBO = F){

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
  fit_s1 = MRAPPSS_EM_func(MRdat,
                          fix.beta = T,
                          beta = 0,
                          pi0 = pi0,
                          sigma.sq = sigma.sq ,
                          tau.sq = tau.sq,
                          Sigma_err = Sigma_err,
                          Omega = Omega,
                          tol = tol,
                          Threshold = Threshold,
                          ELBO = ELBO)

  # stage 2
  fit_s2 = MRAPPSS_EM_func(MRdat,
                         fix.beta = F,
                         beta = 0,
                         pi0 = fit_s1$pi0,
                         sigma.sq = fit_s1$sigma.sq,
                         tau.sq = fit_s1$tau.sq,
                         Sigma_err = Sigma_err,
                         Omega = Omega,
                         tol =  tol,
                         Threshold = Threshold,
                         ELBO = ELBO)
  # Inference
  LR = 2*(fit_s2$likeli-fit_s1$likeli)
  pvalue = pchisq(LR,1,lower.tail = F)
  pvalue = formatC(pvalue, format = "e", digits = 4)
  beta.se = suppressWarnings(abs(fit_s2$beta/sqrt(LR)))
  
  # FBSR
  FBSR = drop(mean(fit_s2$post$Pi * fit_s2$sigma.sq * MRdat$L2)/mean(fit_s2$post$Pi * (Omega[1,1]* MRdat$L2 + Sigma_err[1,1] * MRdat$se.exp^2)))


  cat("***********************************************************\n")
  cat("MR test results of ", exposure , " on ", outcome, ": \n")
  cat("MR-APPSS: beta = ", round(fit_s2$beta,4), "beta.se = ", round(beta.se, 4), "pvalue = ", pvalue, "#SNPs= ", nrow(MRdat), "\n")
  cat("Correlation parameter (rho) due to sample overlap : ", drop(Sigma_err[1,2]), "\n")
  cat("# valid IVs with foreground signals: ", fit_s2$pi0 * nrow(MRdat), "\n")
  cat("Forefround and background signal ratio (FBSR): ", FBSR, "\n")
  cat("***********************************************************\n")

  return( list(MRdat = MRdat,
               exposure = exposure,
               outcome = outcome,
               beta = fit_s2$beta,
               beta.se = beta.se,
               pvalue = pvalue,
               tau.sq = fit_s2$tau.sq,
               sigma.sq = fit_s2$sigma.sq,
               pi0 = fit_s2$pi0,
               post = fit_s2$post,
               FBSR= FBSR,
               likelihoods = fit_s2$likelis,
               Threshold = Threshold,
               method = "MR-APPSS"))
}

