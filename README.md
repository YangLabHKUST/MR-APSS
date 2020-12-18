# MRAPPSS
The MRAPPSS package implement the MR-APPSS approach to test for the causal effects between an exposure and a outcome disease.

The MR-APPSS is a unified approach to Mendelian Randomization accounting for polygenicity, pleiotropy and sample structure using genome-wide summary statistics.
Specifically, MR-APPSS uses a background-foreground model to characterize both SNP-exposure effects and SNP-outcome effects, where the background model accounts for confoundin from genetic correlation and sample structure and the foreground model captures the valid signal for causal inference. 


# Installation 
```{r}
#install.packages("devtools")
devtools::install_github("YangLabHKUST/MR-APSS")
```

# Usage
We illustrate how to analyze GWAS summary level data using the MRAPPSS software by an real example, i.e. BMI (UKB) (exposure) and T2D (outcome). The MRAPSS analysis comprises five steps:

 Step 1: Download GWAS summary-level data from public resources
 
 Step 2: Format data
 
 Step 3: Harmonise datasets and estimate nuisance parameters 
 
 Step 4: IVs selection and LD clumping 
 
 Step 5: Fit MRAPPSS
 
 Step 6: Visualize


The tutorial:  [A real example for perfroming GWAS summary-level data based MR analysis with MRAPPSS package](https://github.com/YangLabHKUST/MRAPSS/blob/master/MRAPSS_Rpackage_Turtorial.pdf) provides details for each step.

To have a quick look at the MRAPPSS, you can skip Steps 1-4 and directly jump to Step 5 and Step 6 to fit MRAPPSS using the outputs we have prepared.
```{r}
library(MRAPPSS)
exposure = "BMI"
outcome = "T2D"
Threshold = 5e-05  # The default p-value threshold IV selection 
data(C)
data(Omega)
data(MRdat)
MRres = MRAPPSS(MRdat,
               exposure="BMI",
               outcome= "T2D",
               C = C,
               Omega =  Omega ,
               Cor.SelectionBias = T)
MRplot(MRres, exposure="BMI", outcome="T2D")
```
The "BMI~T2D" example with 1296 IVs takes about 1 minutes tested on MAC OS 10.14.6 with 1.4 GHz Intel Core i5,16 GB 2133 MHz LPDDR3 and R version 3.6.1. 

# Reproduce[revising]
We provide the downlowad links for 34 GWAS summay datasets used in MR-APPSS paper which can be found [here](https://github.com/YangLabHKUST/MRAPSS_RealData_Code). The R code for processing datasets and implementing other compared MR methods as well as the results can be downloaded [here](https://github.com/YangLabHKUST/MRAPSS_RealData_Code)
# Reference
Xianghong Hu, Jia Zhao, Zhixiang Lin, Yang Wang, Heng Peng, Hongyu Zhao, Xiang Wan, Yang Can, MR-APPSS: a unified approach to Mendelian Randomization accounting for polygenicity, pleiotropy and sample structure using genome-wide summary statistics.

# Development
The MRAPSS package is developed by Xianghong Hu (maxhu@ust.hk).

# Contact information

Please feel free to contact Dr. Xianghong Hu (maxhu@ust.hk) or Prof. Can Yang (macyang@ust.hk) if any questions.
