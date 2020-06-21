# MRAPSS
The MRAPSS package implement the MR-APSS approach to test for the causal effects between an exposure and a outcome disease.

The MR-APSS is a unified approach to Mendelian Randomization accounting for pleiotropy, sample overlap ans selection bias using genome wide summary statistics.

Specifically, MR-APSS uses a background-foreground model to characterize both SNP-exposure effects and SNP-outcome effects, where the background model accounts for the signals due to the shared heritable factors and the foreground model captures the valid signal for causal inference. Building upon the background-foreground model, MR-APSS further takes into account the issues of selection bias and sample overlapping, making it widely applicable for real data analysis.


# Installation 
```{r}
#install.packages("devtools")
devtools::install_github("YangLabHKUST/MRAPSS")
```

# Usage
We illustrate how to analyze GWAS summary level data using the MRAPSS software by an real example, i.e. LDL-C (exposure) and CAD(outcome). The MRAPSS analysis comprises five steps:

 Step 1: Download GWAS summary-level data from public resources
 
 Step 2: Format data
 
 Step 3: Harmonise datasets and estimate nuisance parameters 
 
 Step 4: IVs selection and LD clumping 
 
 Step 5: Fit MRAPSS


The tutorial:  [A real example for perfroming GWAS summary-level data based MR analysis with MRAPSS package](https://github.com/YangLabHKUST/MRAPSS/blob/master/MRAPSS_Rpackage_Turtorial.pdf) provides details for each step.

To have a quick look at the MRAPSS, you can skip Steps 1-4 and directly jump to Step 5 to fit MRAPSS using the outputs we have prepared.
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


# Reference
Xianghong Hu, Jia Zhao, Heng Peng, Yang Wang, Xiang Wan, Yang Can, MR-APSS: a unified approach to Mendelian Randomization accounting for pleiotropy, sample overlap ans selection bias using genome wide summary statistics.

# Development
The MRAPSS package is developed by Xianghong Hu (maxhu@ust.hk).

# Contact information

Please feel free to contact Dr. Xianghong Hu (maxhu@ust.hk) or Prof. Can Yang (macyang@ust.hk) if any questions.
