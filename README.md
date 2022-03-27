# MR-APSS
The MRAPSS package implements the MR-APSS approach to infer the causal relationship between an exposure and an outcome.

MR-APSS is a unified approach to Mendelian Randomization accounting for Pleiotropy and Sample Structure using genome-wide summary statistics. Specifically, MR-APSS uses a foreground-background model to decompose the observed SNP effect sizes, where the background model accounts for confounding factors hidden in GWAS summary statistics, including correlated pleiotropy and sample structure, and the foreground model performs causal inference while accounting for uncorrelated pleiotropy.


# Reproducibility
[Data]
[Real data analysis: negative control outcomes]
[Inferring causal relationships among complex traits]

# Installation 
```{r}
#install.packages("devtools")
devtools::install_github("YangLabHKUST/MR-APSS")
```

# Usage
We illustrate how to analyze GWAS summary level data using the MR-APSS software by a real example, i.e. BMI (UKB) (exposure) and T2D (outcome). The MR-APSS analysis comprises the following steps:

 Step 1: Prepare data and estimate nuisance parameters 
 
 Step 2: Fit MR-APSS for causal inference
 

The tutorial:  [A real example for performing GWAS summary-level data based MR analysis with MRAPSS package](https://github.com/YangLabHKUST/MR-APSS/blob/master/MRAPSS_Rpackage_Tutorial.pdf) provides details for each step.

To have a quick look at the MR-APSS, you can skip Step 1 and directly jump to Step 2 to fit MR-APSS using the outputs we have prepared.
```{r}
library(MRAPSS)
exposure = "BMI"
outcome = "T2D"
Threshold = 5e-05  # The default p-value threshold for IV selection 
data(C)
data(Omega)
data(MRdat)
MRres = MRAPSS(MRdat,
               exposure="BMI",
               outcome= "T2D",
               C = C,
               Omega =  Omega ,
               Cor.SelectionBias = T)
MRplot(MRres, exposure="BMI", outcome="T2D")
```
The "BMI~T2D" example with 1227 IVs takes about 1 minute tested on MAC OS 10.14.6 with 1.4 GHz Intel Core i5,16 GB 2133 MHz LPDDR3 and R version 3.6.1. 

We provide an example R code in "MR-APSS/example" for performing MR analysis with the other five MR methods (IVW, Egger, MRMix, RAPS, and CAUSE). 

# Reference
Xianghong Hu, Jia Zhao, Zhixiang Lin, Yang Wang, Heng Peng, Hongyu Zhao, Xiang Wan, Can Yang. Mendelian Randomization for causal inference accounting for pleiotropy and sample structure using genome-wide summary statistics. bioRxiv 2021.03.11.434915; doi: https://doi.org/10.1101/2021.03.11.434915.

# Contact information

Please feel free to contact Xianghong Hu (maxhu@ust.hk) or Prof. Can Yang (macyang@ust.hk) if any questions.
