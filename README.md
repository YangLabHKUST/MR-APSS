# MR-APSS
<!--[![Github tag](https://badgen.net/github/tag/YangLabHKUST/MR-APSS)](https://github.com/YangLabHKUST/MR-APSS/tags/)
[![Github All Releases](https://img.shields.io/github/downloads/YangLabHKUST/MR-APSS/total.svg)]()
[![Stars](https://img.shields.io/github/stars/YangLabHKUST/MR-APSS?logo=GitHub&color=yellow)](https://github.com/YangLabHKUST/MR-APSS/stargazers) Download counts since July 8, 2022.-->

The MRAPSS package implements the MR-APSS approach to infer the causal relationship between an exposure and an outcome.

MR-APSS is a unified approach to Mendelian Randomization accounting for Pleiotropy and Sample Structure using genome-wide summary statistics. Specifically, MR-APSS uses a foreground-background model to decompose the observed SNP effect sizes, where the background model accounts for confounding factors hidden in GWAS summary statistics, including correlated pleiotropy and sample structure, and the foreground model performs causal inference while accounting for uncorrelated pleiotropy.

![The MR-APSS approach](https://github.com/YangLabHKUST/MRAPSS_RealDataAnalysis_reproduce/blob/master/Fig1.png)
The MR-APSS approach. To infer the causal effect <img src="https://render.githubusercontent.com/render/math?math=\beta"> between exposure X and outcome Y,
MR-APSS uses a foreground-background model to characterize the estimated effects of SNPs <img src="https://render.githubusercontent.com/render/math?math=G_j"> on X and Y (<img src="https://render.githubusercontent.com/render/math?math=\hat\gamma_j"> and <img src="https://render.githubusercontent.com/render/math?math=\hat\Gamma_j">)  with standard errors (<img src="https://render.githubusercontent.com/render/math?math=\hat s_{X,j}"> , <img src="https://render.githubusercontent.com/render/math?math=\hat s_{Y,j}">), where the background model accounts for polygenicity, correlated pleiotropy (**B**) and sample structure (**C**), and the foreground model (**A**) aims to identify informative instruments and account for uncorrelated pleiotropy to perform causal inference. **D** We consider inferring the causal relationship between BMI and T2D as an illustrative example of MR-APSS. The estimated causal effect is indicated by a red line with its 95% confidence interval indicated by the shaded area in transparent red color. Triangles indicate the observed SNP effect sizes (<img src="https://render.githubusercontent.com/render/math?math=\hat\gamma_j"> and <img src="https://render.githubusercontent.com/render/math?math=\hat\Gamma_j">). The color of triangles indicates the posterior of a valid IV, i.e., the posterior of an IV carrying the foreground signal (<img src="https://render.githubusercontent.com/render/math?math=Z_j=1">, dark blue) or not (<img src="https://render.githubusercontent.com/render/math?math=Z_j=0">, light blue)

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

<!-- We provide an example R code in "MR-APSS/example" for performing MR analysis with the other five MR methods (IVW, Egger, MRMix, RAPS, and CAUSE). -->

# Reproducibility
We applied MR-APSS and nine existing summary-level MR methods to (1) test the causal effects of 26 traits on five negative control outcomes (Tanning, Hair color: black, Hair color: blonde; Hair color: dark brown; Hair color: light brown)(130 trait pairs); (2) infer causal relationships between the 26 complex traits (650 trait pairs). In total, there are 780 trait pairs analyzed. We provide [source codes](https://github.com/YangLabHKUST/MRAPSS_RealDataAnalysis_reproduce) for replicating the real data analysis results in the MR-APSS paper. 

**Data download:**  
[Table of GWAS sources](https://github.com/YangLabHKUST/MRAPSS_RealDataAnalysis_reproduce/blob/master/GWAS_26and5_source.csv); [data](https://gohkust-my.sharepoint.com/:u:/g/personal/maxhu_ust_hk/ES9MymhzqMVDhT5L0uwoD0EBHYeLC2wj2CIVcqYm4dJBRQ?e=kvnQlz)(size: 10GB).

**Format data:**  
[code](https://htmlpreview.github.io/?https://github.com/YangLabHKUST/MRAPSS_RealDataAnalysis_reproduce/blob/master/Format_GWASdata.html); [the formatted data](https://gohkust-my.sharepoint.com/:u:/g/personal/maxhu_ust_hk/EVe128zsM5BKm1vze-PaKrUBAaFSvniwNWxkX4HJBO_lJA?e=IRn46T)(size: 711MB).

<!-- We provide [GWAS summary-level datasets] for the five negative control outcomes (Tanning, Hair color: black, Hair color: blonde; Hair color: dark brown; Hair color: light brown) and 26 complex traits. The detailed information for the sources of GWAS datasets is summarized in a [CSV file](https://github.com/YangLabHKUST/MRAPSS_RealDataAnalysis_reproduce/blob/master/GWAS_26and5_source.csv).  
Next, an important step is to format the GWAS summary statistics for the MR_APSS analysis: [code](https://htmlpreview.github.io/?https://github.com/YangLabHKUST/MRAPSS_RealDataAnalysis_reproduce/blob/master/Format_GWASdata.html) and [the formatted datasets](). -->

**Estimate background parameters and plink clumping for the 780 trait pairs:**  
[code for the 130 pairs](https://htmlpreview.github.io/?https://github.com/YangLabHKUST/MRAPSS_RealDataAnalysis_reproduce/blob/master/NC_BackgroundParametersEst_Clumping.html);[code for the 650 pairs](https://htmlpreview.github.io/?https://github.com/YangLabHKUST/MRAPSS_RealDataAnalysis_reproduce/blob/master/Traits_BackgroundParametersEst_Clumping.html); [the estimated background parameters](https://gohkust-my.sharepoint.com/:u:/g/personal/maxhu_ust_hk/ESy5hhtWkpxLt7ikJ55d58QBJxiNqfOUhSEExwVynhPVvA?e=Mf7CjS)(size: 385KB); [data for LD clumped sets of IVs](https://gohkust-my.sharepoint.com/:u:/g/personal/maxhu_ust_hk/EclXsws3J3RNoZDt4Hvgjq4BFPm8amaQuRRpcbt7_eaVHA?e=Agxn2I)(size: 15.1MB).

**Real data analysis with negative control outcomes:**  
[code for MR-APSS](https://htmlpreview.github.io/?https://github.com/YangLabHKUST/MRAPSS_RealDataAnalysis_reproduce/blob/master/NC_MR-APSS.html); [results of MR-APSS](https://github.com/YangLabHKUST/MRAPSS_RealDataAnalysis_reproduce/blob/master/NC_MRAPSS.MRres);  
[code for eight compared methods](https://htmlpreview.github.io/?https://github.com/YangLabHKUST/MRAPSS_RealDataAnalysis_reproduce/blob/master/NC_8MRmethods.html); [results of compared methods](https://github.com/YangLabHKUST/MRAPSS_RealDataAnalysis_reproduce/blob/master/NC_8methods.MRres);  
[code for CAUSE](https://htmlpreview.github.io/?https://github.com/YangLabHKUST/MRAPSS_RealDataAnalysis_reproduce/blob/master/NC_CAUSE.html); [results of CAUSE](https://github.com/YangLabHKUST/MRAPSS_RealDataAnalysis_reproduce/blob/master/NC_CAUSE.MRres);  
[Visualization of results](https://htmlpreview.github.io/?https://github.com/YangLabHKUST/MRAPSS_RealDataAnalysis_reproduce/blob/master/NC-plots.html).
<!-- (https://github.com/YangLabHKUST/MRAPSS_RealDataAnalysis_reproduce/blob/master/NC_CAUSE_MRres). -->
<!-- [code for estimating background parameters and plink clumping](https://htmlpreview.github.io/?https://github.com/YangLabHKUST/MRAPSS_RealDataAnalysis_reproduce/blob/master/NC_BackgroundParametersEst_Clumping.html) and the outputs: [the estimated background parameters](https://gohkust-my.sharepoint.com/:u:/g/personal/maxhu_ust_hk/ESy5hhtWkpxLt7ikJ55d58QBJxiNqfOUhSEExwVynhPVvA?e=Mf7CjS) and [the clumped IV datasets](https://gohkust-my.sharepoint.com/:u:/g/personal/maxhu_ust_hk/ERjc_iN5xbROm3yeZ8fHKtEB6N-0IDcndvbdMrCfztwtsw?e=h5hjpH).  -->

**Inferring causal relationships among complex traits:**  
[code for MR-APSS](https://htmlpreview.github.io/?https://github.com/YangLabHKUST/MRAPSS_RealDataAnalysis_reproduce/blob/master/Traits_MR-APSS.html);[results of MR-APSS](https://github.com/YangLabHKUST/MRAPSS_RealDataAnalysis_reproduce/blob/master/Traits_MRAPSS.MRres);  
[code for eight compared methods](https://htmlpreview.github.io/?https://github.com/YangLabHKUST/MRAPSS_RealDataAnalysis_reproduce/blob/master/Traits_8MRmethods.html); [results of compared methods](https://github.com/YangLabHKUST/MRAPSS_RealDataAnalysis_reproduce/blob/master/Traits_8methods.MRres);  
[code for CAUSE](https://htmlpreview.github.io/?https://github.com/YangLabHKUST/MRAPSS_RealDataAnalysis_reproduce/blob/master/Traits_CAUSE.html); [results of CAUSE](https://github.com/YangLabHKUST/MRAPSS_RealDataAnalysis_reproduce/blob/master/Traits_CAUSE.MRres).   
[Visualization of results](https://htmlpreview.github.io/?https://github.com/YangLabHKUST/MRAPSS_RealDataAnalysis_reproduce/blob/master/Traits-plots.html)
<!-- [code for estimating background parameters and plink clumping](https://htmlpreview.github.io/?https://github.com/YangLabHKUST/MRAPSS_RealDataAnalysis_reproduce/blob/master/Traits_BackgroundParametersEst_Clumping.html) and the outputs: [the estimated background parameters](https://gohkust-my.sharepoint.com/:u:/g/personal/maxhu_ust_hk/ESy5hhtWkpxLt7ikJ55d58QBJxiNqfOUhSEExwVynhPVvA?e=Mf7CjS) and [the clumped IV datasets](https://gohkust-my.sharepoint.com/:u:/g/personal/maxhu_ust_hk/ERjc_iN5xbROm3yeZ8fHKtEB6N-0IDcndvbdMrCfztwtsw?e=h5hjpH);-->

# FAQ
**Q**: What are the quality control criteria for GWAS summary statistics in MR-APSS?  
**A**: MR-APSS uses the following quality control criteria to ensure the quality of data:  
(1). extract SNPs in HapMap 3 list,  
(2). remove SNPs with minor allele frequency < 0.05 (if freq\_col column is available),  
(3). remove SNPs with alleles not in (G, C, T, A),  
(4). remove SNPs with ambiguous alleles (G/C or A/T) or other false alleles (A/A, T/T, G/G or C/C),  
(5). remove SNPs with INFO < 0.9 (if info_col column is available),   
(6). exclude SNPs in the complex Major Histocompatibility Region (Chromosome 6, 26Mb-34Mb),  
(7). remove SNPs with <img src="https://render.githubusercontent.com/render/math?math=\chi^2"> > <img src="https://render.githubusercontent.com/render/math?math=\chi^2_{max}">. The default value for <img src="https://render.githubusercontent.com/render/math?math=\chi^2_{max}"> is max(N/1000, 80).  
 allele frequency in case group and control group
 
**Q**: Does the allele frequency is required for MR-APSS to format GWAS summary statistics? What if I have columns of allele frequency in the case group ("freq_cases") and allele frequency in the control group ("freq_controls")?   
**A**:  Allele frequency ("freq_col"), as well as imputed information ("info_col"), are not required columns in GWAS summary statistics for MR-APSS to obtain the columns in the formatted datasets. To ensure the quality of SNPs, MR-APSS tries to exclude SNPs with  MAF < 0.05.  If allele frequencies are available in summary statistics, MR-APSS will remove SNPs with sample MAF < 0.05.  Considering that sample MAF or imputed information may be missing from the GWAS summary statistics, like LDSC, MR-APSS restricts the analysis to a set of common and well-imputed SNPs in the HapMap 3 reference panel. If both  "freq_cases" and "freq_controls" are available, one can obtain an estimation of MAF using the sample MAF of the control group.  In MR-APSS, one can specify freq_col = “freq_controls”, or ignore "freq_col" when formatting data.

**Q**: What sample size do I need for MR-APSS?  
**A**: In MR-APSS, the background model relies on the assumptions of LDSC (LD Score regression), and the parameters in the background model are estimated using LDSC. Following the practice of [LDSC](https://github.com/bulik/ldsc), we recommend MR-APSS for analyses of GWASs with **more than ~5k samples** to avoid inaccurate estimation results.

**Q**: How does MR-APSS perform LD clumping in real data analysis?  
**A**: In real data analysis, the PLINK LD clumping is used to obtain a subset of nearly independent SNPs as IVs. The default p-value threshold for IV selection for MR-APSS is 5e-05. The squared correlation threshold of clumping (<img src="https://render.githubusercontent.com/render/math?math=r^2">) is 0.001.

**Q**: What exactly does the number "The NO.of valid IVs with foreground signal" reported by MR-APSS mean?  
**A**: The number “NO.of valid IVs with foreground signal” is closely related to the foreground-background model proposed by MR-APSS.  Under the foreground-background model, only a proportion of SNPs with foreground signal (the proportion is denoted by <img src="https://render.githubusercontent.com/render/math?math=\pi_t"> ) will be used for causal inference.  We thus calculated <img src="https://render.githubusercontent.com/render/math?math=\hat\pi_t"> * Total NO. of IVs as “NO.of valid IVs with foreground signal”. This number is also known as the effective number of IVs or the estimated number of valid IVs.

# Reference
<!--Xianghong Hu, Jia Zhao, Zhixiang Lin, Yang Wang, Heng Peng, Hongyu Zhao, Xiang Wan, Can Yang. Mendelian Randomization for causal inference accounting for pleiotropy and sample structure using genome-wide summary statistics. bioRxiv 2021.03.11.434915; doi: https://doi.org/10.1101/2021.03.11.434915. To appear in Proceedings of the National Academy of Sciences (PNAS), 2022.-->
Hu Xianghong, Zhao Jia, Lin Zhixiang, Wang Yang, Peng Heng, Zhao Hongyu, Wan Xiang, and Yang Can. Mendelian randomization for causal inference accounting for pleiotropy and sample structure using genome-wide summary statistics. Proceedings of the National Academy of Sciences, 119(28), July 5, 2022. doi: https://www.pnas.org/doi/10.1073/pnas.2106858119.

# Contact information

Please feel free to contact Xianghong Hu (maxhu@ust.hk) or Prof. Can Yang (macyang@ust.hk) if any questions.
