% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clump.R
\name{clump}
\alias{clump}
\title{Perform LD clumping}
\usage{
clump(
  dat,
  IV.Threshold = 5e-05,
  SNP_col = "SNP",
  pval_col = "pval.exp",
  clump_kb = 1000,
  clump_r2 = 0.001,
  clump_p = 0.999,
  pop = "EUR",
  bfile = NULL,
  plink_bin = NULL
)
}
\arguments{
\item{dat}{a data frame must have columns with information about SNPs and p values}

\item{SNP_col}{column with SNP rsid. The default is `"SNP"`}

\item{pval_col}{column with p value. The default is `"pval"`}

\item{clump_kb}{clumping window in kb. Default is 1000.}

\item{clump_r2}{clumping r2 threshold. Default is 0.001.}

\item{clump_p}{clumping significance level for index variants. Default = 5e-05}

\item{bfile}{bfile as LD reference panel. If this is provided, then will use local PLINK. Default = NULL.}

\item{plink_bin}{path to local plink binary. Default = NULL.}
}
\value{
data frame of clumped SNPs
}
\description{
Peform LD clumping, to prune SNPs in LD within a window. Keep the most significant ones.
}
