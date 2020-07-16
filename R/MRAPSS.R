#'@title  A function for implementing MR-APSS.
#'@description  MR-APSS: a unified approach to Mendelian Randomization accounting for pleiotropy, sample overlap and selection bias  using genome wide summary statistics.
#'MA-APSS uses a variantional EM algorithm for estimation of parameters.
#' MR-APSS uses likelihood ratio test for inference.
#'
#' @param MRdat  data frame at least contain the following varaibles: b.exp b.out se.exp se.out L2. L2:LD score
#' @param exposure exposure name
#' @param outcome   outcome name
#' @param pi0 initial value for pi0, default `NULL` will use the default initialize procedure.
#' @param sigma.sq initial value for sigma.sq , default `NULL`will use the default initialize procedure.
#' @param tau.sq initial value for tau.sq , default `NULL` will use the default initialize procedure.
#' @param Sigma_err the error term correlation matrix. default `diag(2)`.
#' @param Omega  the background varaince component. default `matrix(0,2,2)`.
#' @param tol     tolerence, default '1e-08'
#' @param Threshold   The selection Threshold for correction of selection bias. If Threshold=1, the model won't correct for selection bias.
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
#' \item{method: }{"MR-APSS"}
#' }
#'
#' @examples
#' library(MRAPSS)
#' exposure = "BMI"
#' outcome = "T2D"
#' Threshold = 5e-05  # IV selection Threshold
#' data(Sigma_err)
#' data(Omega)
#' data(MRdat)
#' MRres = MRAPSS(MRdat,
#'                exposure="BMI",
#'                outcome= "T2D",
#'                Sigma_err = Sigma_err,
#'                Omega =  Omega ,
#'                Threshold =  Threshold)
#' MRplot(MRres, exposure="BMI", outcome="T2D")
#' @export

MRAPSS <- function(MRdat=NULL,
                   exposure="exposure",
                   outcome="outcome",
                   pi0=NULL,
                   sigma.sq = NULL,
                   tau.sq = NULL,
                   Sigma_err = matrix(c(1,0,0,1), 2, 2),
                   Omega = matrix(0, 2, 2),
                   tol=1e-08,
                   Threshold=1,
                   ELBO=F){

  if(is.null(MRdat)){
    cat("No data for MR testing")
    return(NULL)
  }

  if(nrow(MRdat) < 4) stop("Not enough IVs.")

  if(is.null(Threshold)) Threshold = max(MRdat$pval.exp)

  if(nrow(Sigma_err)==1 | nrow(Omega) == 1)  jackknife = F

  m = nrow(MRdat)

  # ## stage 0: initialize
  #m_sig = length(which(MRdat$pval.exp <= genome_Threshold))
  #if((is.null(sigma.sq) | is.null(tau.sq)) & m_sig/m > 0.1){

    # initialize pi0, pi0 in c(0.01, 0.99)
    # pi0 = min(m_sig/m, 0.99)
    # pi0 = max(pi0, 0.01)

    # fit_s0 = MRAPSS_EM_func(subset(MRdat, pval.exp<=genome_Threshold),
    #                        pi0 = 1,
    #                        beta =0,
    #                        sigma.sq = sigma.sq,
    #                        tau.sq = tau.sq,
    #                        Sigma_err = Sigma_err,
    #                        Omega = Omega,
    #                        tol = 1e-06,
    #                        Threshold = genome_Threshold)

    # sigma.sq = fit_s0$sigma.sq
    # tau.sq = fit_s0$tau.sq

  #}

  ## stage 1
  fit_s1 = MRAPSS_EM_func(MRdat,
                          fix.beta = T,
                          beta=0,
                          pi0 = pi0,
                          sigma.sq = sigma.sq ,
                          tau.sq = tau.sq,
                          Sigma_err = Sigma_err,
                          Omega = Omega,
                          tol = tol,
                          Threshold = Threshold,
                          ELBO = ELBO)

  # stage 2
  fit_s2 = MRAPSS_EM_func(MRdat,
                         fix.beta = F,
                         beta=0,
                         pi0=fit_s1$pi0,
                         sigma.sq=fit_s1$sigma.sq,
                         tau.sq = fit_s1$tau.sq,
                         Sigma_err = Sigma_err,
                         Omega = Omega,
                         tol =  tol,
                         Threshold = Threshold,
                         ELBO = ELBO)

  LR = 2*(fit_s2$likeli-fit_s1$likeli)
  pvalue = pchisq(LR,1,lower.tail = F)
  pvalue = formatC(pvalue, format = "e", digits = 4)
  beta.se = suppressWarnings(abs(fit_s2$beta/sqrt(LR)))

  cat("***********************************************************\n")
  cat("MR test results of ", exposure , " on ", outcome, ": \n")
  cat("MR-APSS: beta = ", round(fit_s2$beta,4), "beta.se = ", round(beta.se, 4), "pval = ", pvalue, "#SNPs= ", nrow(MRdat), "\n")
  cat("Correlation parameter (rho) due to sample overlap : ", drop(Sigma_err[1,2]), "\n")
  cat("Proportion of effective IVs with forground signals: ", fit_s2$pi0, "\n")
  cat("Variance component (Omega) for background model = \n")
  print(Omega)
  cat("Varaince component (Lambda) for foreground model = \n")
  print(diag(c(fit_s2$sigma.sq, fit_s2$tau.sq)))
  cat("***********************************************************\n")

  return( list(MRdat=MRdat,
               exposure=exposure,
               outcome=outcome,
               beta = fit_s2$beta,
               beta.se = beta.se,
               pvalue = pvalue,
               tau.sq = fit_s2$tau.sq,
               sigma.sq = fit_s2$sigma.sq,
               pi0 = fit_s2$pi0,
               post = fit_s2$post,
               method = "MR-APSS"))
}

