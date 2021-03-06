% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/est_paras.R
\name{est_paras}
\alias{est_paras}
\title{A function harmonising datasets and estimate background parameters by LD score regression.}
\usage{
est_paras(
  dat1,
  dat2,
  trait1.name = "exposure",
  trait2.name = "outcome",
  LDSC = T,
  h2.fix.intercept = F,
  ldscore.dir = NULLL
)
}
\arguments{
\item{dat1:}{formmated summary statistics for trait 1.}

\item{dat2:}{formmated summary statistics for trait 2.}

\item{trait1.name:}{specify the name of trait 1, default `exposure`.}

\item{trait2.name:}{specify the name of trait 2, default `outcome`.}

\item{LDSC:}{whether to run LD score regression, default `TRUE`. If `FALSE`, the function will not give the parameter estimates but will do harmonising.}

\item{h2.fix.intercept:}{whether to fix LD score regression intercept to 1, default `FALSE`.}

\item{ldscore.dir:}{specify the path to the LD score files.}
}
\value{
List with the following elements:
\describe{
\item{Mdat}{Homonised data set }
\item{C}{the estimated C matrix capturing the effects of sample structure }
\item{Omega}{the estimated variance-covariance matrix for polygenic effects}
}
}
\description{
A function harmonising datasets and estimate background parameters by LD score regression.
}
