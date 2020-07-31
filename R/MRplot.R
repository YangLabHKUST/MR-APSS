#' @title  Visualize the MRAPSS results
#'
#' @param MRres: MRAPSS fit results
#' @param exposure: exposure name
#' @param outcome : outcome name
#'
#' @return Plot of SNP-exposure effect and SNP-outcome effect with the causal effect and 95% confidence interval.
#' @export
#'
MRplot <- function(MRres, exposure="trait 1", outcome ="trait 2"){
  MRdat = MRres$MRdat

  x =ifelse(abs(MRdat$b.exp - MRdat$se.exp) > abs(MRdat$b.exp + MRdat$se.exp),
            MRdat$b.exp - MRdat$se.exp-0.01, MRdat$b.exp + MRdat$se.exp+0.01)


  library(ggnewscale)
  ggplot2::ggplot(MRdat, ggplot2::aes(x=b.exp, y=b.out)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = b.out - se.out, ymax = b.out + se.out),colour="gray60", size=0.2) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = b.exp - se.exp, xmax = b.exp + se.exp, height = 0),colour="gray60", size=0.2
    ) +
    ggplot2::geom_point(data=MRdat, ggplot2::aes(x=b.exp, y=b.out, size = -log10(pval.exp), color=MRres$post$Pi), shape=17) +
    ggplot2::scale_color_continuous(name = "IV Posterior", high = "#132B43", low = "#56B1F7")+
    ggplot2::geom_abline(ggplot2::aes(slope = MRres$beta, intercept=0), color= "red", size=0.5,alpha=1) +
    ggplot2::geom_hline(yintercept = 0, lty="dotted", color="gray60") +
    ggplot2::labs(x= paste("SNP effect on", exposure), y= paste("SNP effect on", outcome)) +
    ggplot2::geom_ribbon(ggplot2::aes(x = x,
                    ymin = (MRres$beta -1.96*MRres$beta.se)*x,
                    ymax = (MRres$beta +1.96*MRres$beta.se)*x), fill="red",alpha=0.15) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.border = ggplot2::element_blank(),
          panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          # legend.direction = "horizontal",
          # legend.justification=c(0,1),
          # legend.position=c(0.015,1),
          legend.position = "right",
          legend.text = ggplot2::element_text(hjust = 0, size=15),
          legend.title = ggplot2::element_text(hjust = 0, size=15),
          #legend.title = element_blank(),
          axis.text = ggplot2::element_text(size=10, color = "black"),
          axis.title = ggplot2::element_text(size=15, face = "bold", color = "black"),
          plot.title = ggplot2::element_text(hjust=-0.15, size = 20, face = "bold"),
          plot.subtitle = ggplot2::element_text(vjust=0, hjust=0.5, size = 15, face = "bold"),
          axis.line = ggplot2::element_line(colour = "black"))


}
