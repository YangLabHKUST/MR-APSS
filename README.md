# MRAPSS
The MRAPSS package implement the MR-APSS approach to test for the causal effects between an exposure and a outcome disease.

The MR-APSS is a unified approach to Mendelian Randomization accounting for pleiotropy, sample overlap ans selection bias using genome wide summary statistics. 

This R package is developed by HU Xianghong(maxhu@ust.hk).

# Reference
Xianghong Hu, Jia Zhao, Heng Peng, Yang Wang, Xiang Wan, Yang Can, MR-APSS: a unified approach to Mendelian Randomization accounting for pleiotropy, sample overlap ans selection bias using genome wide summary statistics.

# Installation 
```{r}
#install.packages("devtools")
devtools::install_github("YangLabHKUST/MRAPSS")
```

# An example
To have a quick look at the MRAPSS, run 
```{r}
library(MRAPSS)
exposure = "LDL-C"
outcome = "CAD"
Threshold = 5e-05  # IV selection Threshold
data(Omega)        # estimates of covariance matrix of background effects
data(Sigma_err)    # estimates of correlation matrix for error term
data(MRdat)        # dataset for clumped SNPs
MRres = MRAPSS(MRdat,
               exposure = "LDl-C",
               outcome = "CAD",
               Sigma_err = Sigma_err,
               Omega =  Omega ,
               Threshold =  Threshold)
MRplot(MRres, exposure="LDL-C", outcome="CAD")
```

See [A real example for perfroming GWAS summary-level data based MR analysis with MRAPSS package](https://github.com/hxh0504/MRAPSS/blob/master/Turtorial.pdf) for details.


