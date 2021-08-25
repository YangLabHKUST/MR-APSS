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


  cat("***********************************************************\n")
  cat("MR test results of ", exposure , " on ", outcome, ": \n")
  cat("r=", r, "beta = ", round(fit_s2$beta,4), "beta.se = ", round(beta.se, 4), "pvalue = ", pvalue, "#SNPs= ", nrow(MRdat), "\n")
  cat("# valid IVs with foreground signals: ", fit_s2$pi0 * nrow(MRdat), "\n")
  cat("***********************************************************\n")
  return( list(exposure = exposure,
               outcome = outcome,
               beta = fit_s2$beta,
               beta.se = beta.se,
               pvalue = pvalue,
               Threshold = Threshold))
}


