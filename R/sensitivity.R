#' @export
sensitivity <- function(MRdat=NULL,
                        Omega = NULL,
                        C=NULL,
                        exposure="exposure",
                        outcome ="outcome"){

  fit0 = MR_EM_fix_rf_func(MRdat,
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
    fit = MR_fix_rf(MRdat,
                  r=r,
                  C = C,
                  Omega =  Omega ,
                  Cor.SelectionBias = T,
                  tol=1e-12)

    res = rbind(res, data.frame(r=r, beta= fit$beta, se=fit$beta.se, pval = fit$pvalue))
  }
    res$r = round(res$r, 3)
    plot = ggplot2::ggplot(data=res) +
      ggplot2::geom_point(ggplot2::aes(x= r, y=beta), shape=16,size=5) +
      ggplot2::geom_errorbar(ggplot2::aes(x=r, y=beta, ymin = beta - 1.96*se, ymax = beta + 1.96*se), width=0.01) +
      ggplot2::geom_line(ggplot2::aes(x= r, y=beta, group=1), lty="dotted") +
      ggplot2::geom_hline(yintercept = 0, lty = "dotted") +
      ggplot2::labs(x = "Correlation between IV strength and direct effect in the foreground model",
           y = "Causal effect estimate",
           title = paste0( exposure, " and ", outcome)) +
      ggplot2::scale_x_continuous(breaks = res$r)+
      ggplot2::theme_classic() +
      ggplot2::theme(axis.text = ggplot2::element_text(size= 12),
                   axis.title = ggplot2::element_text(size= 12),
            plot.title = ggplot2::element_text(size= 25))

    return(list(estimates = res,
                plot = plot))

}


