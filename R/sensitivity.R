#' @export
sensitivity <- function(MRdat=NULL,
                        Omega = NULL,
                        C=NULL){

  fit0 = MR_EM_fixr_func(MRdat,
                         fix.beta = T,
                         beta = 0,
                         r=0,
                         fix.r = F,
                         pi0 = NULL,
                         Sigma=NULL,
                         C = C,
                         Omega = Omega,
                         tol = 1e-08,
                         Threshold =  unique(MRdat$Threshold))

  rmax = floor(fit0$Sigma[1,2]/sqrt(fit0$Sigma[1,1] * fit0$Sigma[2,2])*100)/100

  res= NULL
  for( r in seq(0, rmax, length.out = 10)){
    fit = MR_fixr(MRdat,
                  r=r,
                  C = C,
                  Omega =  Omega ,
                  Cor.SelectionBias = T,
                  tol=1e-12)

    res = rbind(res, data.frame(r=r, beta= fit$beta, se=fit$beta.se, pval = fit$pvalue))
  }

    plot = ggplot(data=res) +
      geom_point(aes(x= r, y=beta), shape=16,size=5) +
      geom_errorbar(aes(x=r, y=beta, ymin = beta - 1.96*se, ymax = beta + 1.96*se), width=0.01) +
      geom_line(aes(x= r, y=beta, group=1), lty="dotted") +
      geom_hline(yintercept = 0, lty = "dotted") +
      labs(x = TeX("$r_f$ "),
           y = "Causal effect estimate",
           title = paste0( exp.name, " and ", out.name)) +
      scale_x_continuous(breaks = res$r)+
      theme_classic() +
      theme(axis.text = element_text(size= 12),
            axis.title = element_text(size= 20),
            plot.title = element_text(size= 25))

    return(list(estimates = res,
                plot = plot))

}


