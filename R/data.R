# # w_hm3.snplist: the hapmap3 snplist downloaded from https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2
# w_hm3.snplist = readr::read_delim("data-raw/w_hm3.snplist", '\t', escape_double = FALSE,
#                         trim_ws = TRUE, progress = F, col_names =T)
#
# # MHC.SNPs : SNPs in MHC region
# MHC.SNPs = readr::read_delim("data-raw/MHC.SNPs", '\t', escape_double = FALSE,
#                          trim_ws = TRUE, progress = F, col_names =T)$SNP
#
# # # Sigma_err estimate for LDL-C~CAD example
# Sigma_err = read.table("data-raw/LDL-C~CAD_Sigma_err", header = F)
# Sigma_err = as.matrix(Sigma_err,2,2)
#
# # # Omega estimate for LDL-C~CAD example
# Omega = read.table("data-raw/LDL-C~CAD_Omega", header = F)
# Omega = as.matrix(Omega,2,2)
#
# # # Dataset for MRAPSS analysis (LDL-C~CAD example)
# MRdat = read.table("data-raw/LDL-C~CAD", header = T)
#
# # # Dataset for MRAPSS analysis (CAD~LDL-C example)
# MRdat_rev = read.table("data-raw/CAD~LDL-C", header = T)
#
# usethis::use_data(w_hm3.snplist, MHC.SNPs, Sigma_err, Omega, MRdat, MRdat_rev, overwrite = TRUE)

