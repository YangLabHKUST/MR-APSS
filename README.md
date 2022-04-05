# MR-APSS
The MRAPSS package implements the MR-APSS approach to infer the causal relationship between an exposure and an outcome.

MR-APSS is a unified approach to Mendelian Randomization accounting for Pleiotropy and Sample Structure using genome-wide summary statistics. Specifically, MR-APSS uses a foreground-background model to decompose the observed SNP effect sizes, where the background model accounts for confounding factors hidden in GWAS summary statistics, including correlated pleiotropy and sample structure, and the foreground model performs causal inference while accounting for uncorrelated pleiotropy.

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
The "BMI~T2D" example with 1227 IVs takes about 1 minute when tested on MAC OS 10.14.6 with 1.4 GHz Intel Core i5,16 GB 2133 MHz LPDDR3 and R version 3.6.1. 

We provide an example R code in "MR-APSS/example" for performing MR analysis with the other five MR methods (IVW, Egger, MRMix, RAPS, and CAUSE). 

# Reproducibility (in progress)
We applied MR-APSS and nine existing summry-level MR methods to: (1) test the causal effects of 26 traits on five negative control outcomes (Tanning, Hair color: black, Hair color: blonde; Hair color: dark brown; Hair color: light brown)(130 trait pairs); (2) infer causal relationships between the 26 complex traits (650 trait pairs). In total, there are 780 trait pairs analyzed. We provide the [sources codes](https://github.com/YangLabHKUST/MRAPSS_RealDataAnalysis_reproduce) for replicating the real data analysis results in the MR-APSS paper. 

**Data download:**  
[Table of GWAS sources](https://github.com/YangLabHKUST/MRAPSS_RealDataAnalysis_reproduce/blob/master/GWAS_26and5_source.csv); [data]().

**Format data:**  
[code](https://htmlpreview.github.io/?https://github.com/YangLabHKUST/MRAPSS_RealDataAnalysis_reproduce/blob/master/Format_GWASdata.html); [the formatted data](https://gohkust-my.sharepoint.com/:u:/g/personal/maxhu_ust_hk/EVe128zsM5BKm1vze-PaKrUBAaFSvniwNWxkX4HJBO_lJA?e=2C874o).

<!-- We provide [GWAS summary-level datasets] for the five negative control outcomes (Tanning, Hair color: black, Hair color: blonde; Hair color: dark brown; Hair color: light brown) and 26 complex traits. The detailed information for the sources of GWAS datasets is summarized in a [CSV file](https://github.com/YangLabHKUST/MRAPSS_RealDataAnalysis_reproduce/blob/master/GWAS_26and5_source.csv).  
Next, a important step is to format the GWAS summay statistics for using in the MR_APSS analysis: [code](https://htmlpreview.github.io/?https://github.com/YangLabHKUST/MRAPSS_RealDataAnalysis_reproduce/blob/master/Format_GWASdata.html) and [the formatted datasets](). -->

**Estimate background parameters and plink clumping for the 780 trait pairs:**  
[code](https://htmlpreview.github.io/?https://github.com/YangLabHKUST/MRAPSS_RealDataAnalysis_reproduce/blob/master/NC_BackgroundParametersEst_Clumping.html); [the estimated background parameters](https://gohkust-my.sharepoint.com/:u:/g/personal/maxhu_ust_hk/ESy5hhtWkpxLt7ikJ55d58QBJxiNqfOUhSEExwVynhPVvA?e=Mf7CjS); [the clumped IV datasets](https://gohkust-my.sharepoint.com/:u:/g/personal/maxhu_ust_hk/ERjc_iN5xbROm3yeZ8fHKtEB6N-0IDcndvbdMrCfztwtsw?e=h5hjpH).  

**Real data analysis with negative control outcomes:**  
[code for MR-APSS](https://htmlpreview.github.io/?https://github.com/YangLabHKUST/MRAPSS_RealDataAnalysis_reproduce/blob/master/NC_MR-APSS.html) and [results of MR-APSS](https://github.com/YangLabHKUST/MRAPSS_RealDataAnalysis_reproduce/blob/master/NC_MRAPSS.MRres);  
[code for eight compared methods] and [results of compared methods];  
[code for CAUSE](https://htmlpreview.github.io/?https://github.com/YangLabHKUST/MRAPSS_RealDataAnalysis_reproduce/blob/master/NC_CAUSE.html) and [results of CAUSE](https://github.com/YangLabHKUST/MRAPSS_RealDataAnalysis_reproduce/blob/master/NC_CAUSE_MRres). 
<!-- [code for estimating background parameters and plink clumping](https://htmlpreview.github.io/?https://github.com/YangLabHKUST/MRAPSS_RealDataAnalysis_reproduce/blob/master/NC_BackgroundParametersEst_Clumping.html) and the outputs: [the estimated background parameters](https://gohkust-my.sharepoint.com/:u:/g/personal/maxhu_ust_hk/ESy5hhtWkpxLt7ikJ55d58QBJxiNqfOUhSEExwVynhPVvA?e=Mf7CjS) and [the clumped IV datasets](https://gohkust-my.sharepoint.com/:u:/g/personal/maxhu_ust_hk/ERjc_iN5xbROm3yeZ8fHKtEB6N-0IDcndvbdMrCfztwtsw?e=h5hjpH).  -->

**Inferring causal relationships among complex traits:**  
[code for MR-APSS](https://htmlpreview.github.io/?https://github.com/YangLabHKUST/MRAPSS_RealDataAnalysis_reproduce/blob/master/Traits_MR-APSS.html) and [results of MR-APSS](https://github.com/YangLabHKUST/MRAPSS_RealDataAnalysis_reproduce/blob/master/Traits_MRAPSS.MRres);  
[code for eight compared methods] and [results of compared methods];  
[code for CAUSE] and [results of CAUSE]. 

<!-- [code for estimating background parameters and plink clumping](https://htmlpreview.github.io/?https://github.com/YangLabHKUST/MRAPSS_RealDataAnalysis_reproduce/blob/master/Traits_BackgroundParametersEst_Clumping.html) and the outputs: [the estimated background parameters](https://gohkust-my.sharepoint.com/:u:/g/personal/maxhu_ust_hk/ESy5hhtWkpxLt7ikJ55d58QBJxiNqfOUhSEExwVynhPVvA?e=Mf7CjS) and [the clumped IV datasets](https://gohkust-my.sharepoint.com/:u:/g/personal/maxhu_ust_hk/ERjc_iN5xbROm3yeZ8fHKtEB6N-0IDcndvbdMrCfztwtsw?e=h5hjpH);-->


# Reference
Xianghong Hu, Jia Zhao, Zhixiang Lin, Yang Wang, Heng Peng, Hongyu Zhao, Xiang Wan, Can Yang. Mendelian Randomization for causal inference accounting for pleiotropy and sample structure using genome-wide summary statistics. bioRxiv 2021.03.11.434915; doi: https://doi.org/10.1101/2021.03.11.434915.

# Contact information

Please feel free to contact Xianghong Hu (maxhu@ust.hk) or Prof. Can Yang (macyang@ust.hk) if any questions.
