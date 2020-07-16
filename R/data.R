# # w_hm3.snplist: the hapmap3 snplist downloaded from https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2
# w_hm3.snplist = readr::read_delim("data-raw/w_hm3.snplist", '\t', escape_double = FALSE,
#                         trim_ws = TRUE, progress = F, col_names =T)
#
# # MHC.SNPs : SNPs in MHC region
# MHC.SNPs = readr::read_delim("data-raw/MHC.SNPs", '\t', escape_double = FALSE,
#                          trim_ws = TRUE, progress = F, col_names =T)$SNP
#
# # # Sigma_err estimate for BMI~T2D example
# Sigma_err = read.table("data-raw/BMI~T2D_Sigma_err", header = F)
# Sigma_err = as.matrix(Sigma_err,2,2)
#
# # # Omega estimate for BMI~T2D example
# Omega = read.table("data-raw/BMI~T2D_Omega", header = F)
# Omega = as.matrix(Omega,2,2)
#
# # # Dataset for MRAPSS analysis (BMI~T2D example)
# MRdat = read.table("data-raw/BMI~T2D", header = T)
#
# # # Dataset for MRAPSS analysis (T2D ~ BMI example)
# ##MRdat_rev = read.table("data-raw/t2D ~ BMI", header = T)
#
# usethis::use_data(w_hm3.snplist, MHC.SNPs, Sigma_err, Omega, MRdat, overwrite = TRUE)
#
