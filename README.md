# MRAPSS
The MRAPSS package implement the MR-APSS approach to test for the causal effects between an exposure and a outcome disease.

The MR-APSS is a unified approach to Mendelian Randomization accounting for pleiotropy, sample overlap ans selection bias using genome wide summary statistics.

Specifically, MR-APSS uses a background-foreground model to characterize both SNP-exposure effects and SNP-outcome effects, where the background model accounts for the signals due to the shared heritable factors and the foreground model captures the valid signal for causal inference. Building upon the background-foreground model, MR-APSS further takes into account the issues of selection bias and sample overlapping, making it widely applicable for real data analysis.


# Installation 
```{r}
#install.packages("devtools")
devtools::install_github("YangLabHKUST/MRAPSS")
```

# An example
To have a quick look at the MRAPSS, run the following example using the data we have prepared.
```{r}
library(MRAPSS)
exposure = "LDL-C"
outcome = "CAD"
Threshold = 5e-05  # IV selection Threshold
data(Omega)        # estimates of covariance matrix of background effects
data(Sigma_err)    # estimates of correlation matrix for error term
data(MRdat)        # dataset for clumped SNPs
MRres = MRAPSS(MRdat,
               exposure = "LDL-C",
               outcome = "CAD",
               Sigma_err = Sigma_err,
               Omega =  Omega ,
               Threshold =  Threshold)
MRplot(MRres, exposure="LDL-C", outcome="CAD")
```

To apply MR-APSS with GWAS summary-level data, please see  [A real example for perfroming GWAS summary-level data based MR analysis with MRAPSS package](https://github.com/hxh0504/MRAPSS/blob/master/Turtorial.pdf) for details.

# Reference
Xianghong Hu, Jia Zhao, Heng Peng, Yang Wang, Xiang Wan, Yang Can, MR-APSS: a unified approach to Mendelian Randomization accounting for pleiotropy, sample overlap ans selection bias using genome wide summary statistics.

# Developer
This R package is developed and maintained by HU Xianghong(maxhu@ust.hk). 
