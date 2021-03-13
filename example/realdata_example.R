install.packages("devtools")
library(devtools)
install_github("MRCIEU/TwoSampleMR")
install_github("qingyuanzhao/mr.raps")
devtools::install_github("gqi/MRMix")
devtools::install_github("jean997/cause")

library(MRAPSS)
exposure = "BMI"
outcome = "T2D"
Threshold = 5e-05  # IV selection Threshold
data(C)
data(Omega)
data(MRdat)

# MR-APSS
MRres = MRAPSS(MRdat,
               exposure = "BMI",
               outcome = "T2D",
               C = C,
               Omega =  Omega ,
               Cor.SelectionBias = T)
MRplot(MRres, exposure = "BMI", outcome = "T2D")

# MR-APSS(Omega=0)
MRres = MRAPSS(MRdat,
               exposure = "BMI",
               outcome = "T2D",
               C = C,
               Cor.SelectionBias = T)

# MR-APSS(C=I)
MRres = MRAPSS(MRdat,
               exposure = "BMI",
               outcome = "T2D",
               Omega =  Omega ,
               Cor.SelectionBias = T)

# IVW
indx =  which(MRdat$pval.exp <=  5e-08)
res.IVW <- TwoSampleMR::mr_ivw_fe(MRdat$b.exp[indx], MRdat$b.out[indx],
                                  MRdat$se.exp[indx], MRdat$se.out[indx])


# RAPS
df = data.frame(beta.exposure = MRdat$b.exp[indx], 
                beta.outcome = MRdat$b.out[indx],
                se.exposure = MRdat$se.exp[indx],
                se.outcome = MRdat$se.out[indx])
res.raps <- mr.raps::mr.raps(df, diagnostics=F)

# MRMix 
res.MRMix = MRMix::MRMix(MRdat[indx,]$b.exp, MRdat[indx,]$b.out, MRdat[indx,]$se.exp, MRdat[indx,]$se.out)

# Egger
res.egger <- TwoSampleMR::mr_egger_regression(MRdat$b.exp[indx], MRdat$b.out[indx],
                                              MRdat$se.exp[indx], MRdat$se.out[indx])

# CAUSE
# You should first download the Summary statistics for BMI and T2D at the following link.
# https://gohkust-my.sharepoint.com/:f:/g/personal/maxhu_ust_hk/Egqda8NZeKdEu-7lg2iPqOMBKTskTogtjeysUW9y8TXrfA?e=hqHbvw

Threshold=1e-03
# read in the datasets
X1 = readr::read_delim(BMI.dir,delim="\t",
                       escape_double = FALSE,
                       trim_ws = TRUE,
                       progress = F)

X2 = readr::read_delim(T2D.dir, "\t",
                       escape_double = FALSE,
                       trim_ws = TRUE,
                       progress = F)

X1$b = X1$Z/sqrt(X1$N)
X2$b = X2$Z/sqrt(X2$N)
X1$se = 1/sqrt(X1$N)
X2$se = 1/sqrt(X2$N)

X <- try(cause::gwas_merge(X1, X2, 
                    snp_name_cols = c("SNP", "SNP"),
                    beta_hat_cols = c("b", "b"),
                    se_cols = c("se", "se"),
                    A1_cols = c("A1", "A1"),
                    A2_cols = c("A2", "A2")))

d0 = X1[, c("SNP", "P")]
colnames(d0) = c("snp", "pval.exp")
X0 = merge(X, d0, by="snp")

# The clumped dataset "clumped_3" is avaliable in ./example/BMI_ukb~T2D_CAUSE_clumped.RData
clumped_3 = MRAPSS::clump(X0,
                          IV.Threshold = 1e-03,
                          SNP_col = "snp",
                          pval_col = "pval.exp",
                          clump_kb = 1000,
                          clump_r2 = 0.1)

varlist <- with(X, sample(snp, size=min(nrow(X), 1000000), replace=FALSE))

# "params" is avaliable in ./example/BMI_ukb~T2D_CAUSE_paras.RData
params <- try(cause::est_cause_params(X, varlist))

top_ldl_pruned_vars =intersect(as.character(X$snp), as.character(subset(clumped, pval.exp <= Threshold)$snp))
  
cause_res <- try(cause::cause(X=X, variants = top_ldl_pruned_vars , param_ests = params, force=TRUE))
  
res_elpd <- data.frame(trait1.name,
                         trait2.name,
                         Threshold,
                         length(top_ldl_pruned_vars),
                         cause_res$elpd)
res.cause.est = summary(cause_res, ci_size=0.95)
  

