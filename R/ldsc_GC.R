ldsc_GC <- function(merged,
                    trait1.name = "exposure",
                    trait2.name = "outcome",
                    M,
                    Twostep = T,
                    h2.fix.intercept = T,
                    n.blocks = 200){

  #time <- proc.time()

  # Storage:
  cov <- matrix(NA,nrow=2,ncol=2)
  cov.se <- matrix(NA,nrow=2,ncol=2)
  I <- matrix(NA,nrow=2,ncol=2)
  I.se <- matrix(NA,nrow=2,ncol=2)

  merged <- merged[with(merged,order(CHR,BP)),]
  n.snps = nrow(merged)

  merged$Zxy = merged$Z.x*merged$Z.y
  merged <- merged[,c("SNP","chi2.x","chi2.y","N.x","N.y","Zxy","L2")]

  merged11 = merged[,c("SNP","chi2.x","chi2.x","N.x","N.x","chi2.x","L2")]
  colnames(merged11) = c("SNP","chi2.x","chi2.y","N.x","N.y","Zxy","L2")

  merged22 = merged[,c("SNP","chi2.y","chi2.y","N.y","N.y","chi2.y","L2")]
  colnames(merged22) = c("SNP","chi2.x","chi2.y","N.x","N.y","Zxy","L2")

  if(Twostep){
    step1.idx = which(merged$chi2.x<30 & merged$chi2.y <30)
    seperator <- floor(seq(from=1, to=length(step1.idx), length.out =(n.blocks+1)))
    new.seperator <-  c(1, step1.idx[seperator[2:n.blocks]], nrow(merged))
  }else{
    step1.idx = NULL
    new.seperator <- floor(seq(from=1, to=nrow(merged), length.out =(n.blocks+1)))
  }

  # estmate heritability
  message("Estimate heritability for trait 1 ...")

  if(h2.fix.intercept){
    h2.1 = est_h2(merged = merged[, c("SNP","chi2.x","N.x","L2")], M, trait1.name, Twostep, h2.fix.intercept, n.blocks=n.blocks, step1.idx, jknife=F)
    h2.2 = est_h2(merged = merged[, c("SNP","chi2.y","N.y","L2")], M, trait2.name, Twostep, h2.fix.intercept, n.blocks=n.blocks, step1.idx, jknife=F)
    h2.1F = est_h2(merged = merged[, c("SNP","chi2.x","N.x","L2")], M, trait1.name, Twostep, F, step1.idx,n.blocks=n.blocks, jknife=F)
    h2.2F = est_h2(merged = merged[, c("SNP","chi2.y","N.y","L2")], M, trait2.name, Twostep, F,  step1.idx, n.blocks=n.blocks, jknife=F)
  }else{
    h2.1 = est_h2(merged = merged[, c("SNP","chi2.x","N.x","L2")], M, trait1.name, Twostep, h2.fix.intercept, step1.idx, n.blocks=n.blocks, jknife=F)
    h2.2 = est_h2(merged = merged[, c("SNP","chi2.y","N.y","L2")], M, trait2.name, Twostep, h2.fix.intercept, step1.idx, n.blocks=n.blocks, jknife=F)
    h2.1F = h2.1
    h2.2F = h2.2
  }

  if(h2.fix.intercept){
    h2_1 = est_gencov(merged = merged11,
                      h1 = h2.1$h2,
                      h2 = h2.1$h2,
                      intercept.h1 = h2.1$intercept,
                      intercept.h2 = h2.1$intercept,
                      intercept.gencov = h2.1$intercept,
                      n.blocks=n.blocks,
                      M,
                      Twostep,
                      fix.gcov.intercept=T,
                      step1.idx,jknife=T)
  }else{
    h2_1 = est_gencov(merged = merged11,
                      h1 = h2.1F$h2,
                      h2 = h2.1F$h2,
                      intercept.h1 = h2.1F$intercept,
                      intercept.h2 = h2.1F$intercept,
                      n.blocks=n.blocks,
                      intercept.gencov = h2.1F$intercept,
                      M,Twostep,
                      fix.gcov.intercept=F,
                      step1.idx,jknife=T)
  }

  message("Mean Chi2:",round(mean(merged$chi2.x),4),".")
  message("Intercept: ",round(h2_1$intercept,4),"(",round(h2_1$intercept.se,4),").")
  message("Total Observed scale h2:",round(h2_1$rho_g,4),"(",round(h2_1$rho_g.se,4),").\n")

  message("Estimate heritability for trait 2 ...")
  if(h2.fix.intercept){
    h2_2 = est_gencov(merged = merged22,
                      h1 = h2.2$h2,
                      h2 = h2.2$h2,
                      intercept.h1 = h2.2$intercept,
                      intercept.h2 = h2.2$intercept,
                      n.blocks=n.blocks,
                      intercept.gencov = h2.2$intercept,
                      M,Twostep,fix.gcov.intercept=T,
                      step1.idx,jknife=T)
  }else{
    h2_2 = est_gencov(merged = merged22,
                      h1 = h2.2F$h2,
                      h2 = h2.2F$h2,
                      intercept.h1 = h2.2F$intercept,
                      intercept.h2 = h2.2F$intercept,
                      n.blocks=n.blocks,
                      intercept.gencov = h2.2F$intercept,
                      M,Twostep,fix.gcov.intercept=F,
                      step1.idx,jknife=T)
  }

  message("Mean Chi2:",round(mean(merged$chi2.y),4), ".")
  message("Intercept: ",round(h2_2$intercept,4)," (",round(h2_2$intercept.se,4),").")
  message("Total Observed scale h2:",round(h2_2$rho_g,4)," (",round(h2_2$rho_g.se,4),").\n")


  # estimate genetic covariance
  message("Estimate genetic covariance ...")
  if(h2.fix.intercept){
    rho_g = est_gencov(merged = merged,
                       h1 = h2.1$h2,
                       h2 = h2.2$h2,
                       intercept.h1 = h2.1$intercept,
                       intercept.h2 = h2.2$intercept,
                       n.blocks=n.blocks,
                       intercept.gencov=0,
                       M,Twostep,fix.gcov.intercept=F,
                       step1.idx,jknife=T)
  }else{
    rho_g = est_gencov(merged = merged,
                       h1 = h2.1F$h2,
                       h2 = h2.2F$h2,
                       intercept.h1 = h2.1F$intercept,
                       intercept.h2 = h2.2F$intercept,
                       n.blocks=n.blocks,
                       intercept.gencov=0,
                       M,Twostep,fix.gcov.intercept=F,
                       step1.idx,jknife=T)
  }

  message("Intercept: ",round(rho_g$intercept,4)," (",round(rho_g$intercept.se,4),").")
  message("Total Observed scale gencov: ",round(rho_g$rho_g,4)," (",round(rho_g$rho_g.se,4),").")


  cov[1,1] = h2_1$rho_g
  cov[2,2] = h2_2$rho_g
  I[1,1] = h2_1$intercept
  I[2,2] = h2_2$intercept

  cov[1,2] <- cov[2,1] <- rho_g$rho_g
  I[1,2] <-  I[2,1] <- rho_g$intercept
  rg = rho_g$rho_g/sqrt(h2_1$rho_g * h2_2$rho_g)


  cov.se[1,1] <- h2_1$rho_g.se
  cov.se[2,2] <- h2_2$rho_g.se
  cov.se[1,2]<- cov.se[2,1]<- rho_g$rho_g.se

  I.se[1,1] <- h2_1$intercept.se
  I.se[2,2] <- h2_2$intercept.se
  I.se[1,2] <-I.se[2,1] <- rho_g$intercept.se

  #time_all <- proc.time()-time

  #message("Time elapsed: ", time_all[3])

  denome_delete.values <- sqrt(h2_1$delete.values[,1] * h2_2$delete.values[,1])

  pseudo.values = n.blocks*rg - (n.blocks-1)* rho_g$delete.values[,1]/denome_delete.values

  jackknife.cov <- var(pseudo.values)/n.blocks

  rg.se <- sqrt(jackknife.cov)

  # if(h2.fix.intercept){
  #   h2_1_inter_delete.values = 1
  #   h2_2_inter_delete.values = 1
  # }else{
  #   h2_1_inter_delete.values = h2_1$delete.values[,2]
  #   h2_2_inter_delete.values = h2_2$delete.values[,2]
  # }
  #
  # h2_1_delete.values = h2_1$delete.values[,1]
  # h2_2_delete.values = h2_2$delete.values[,1]
  # gcov_delete.values = rho_g$delete.values[,1]
  # rho_delete.values = rho_g$delete.values[,2]
  #
  # Sigma_LD_delete = as.matrix(data.frame(h2_1_inter_delete.values,
  #                                        rho_delete.values,
  #                                        rho_delete.values,
  #                                        h2_2_inter_delete.values))
  #
  # Omega_delete = as.matrix(data.frame(h2_1_delete.values,
  #                                     gcov_delete.values,
  #                                     gcov_delete.values,
  #                                     h2_2_delete.values))

  return(list(cov = cov,
              cov.se = cov.se,
              I = I,
              I.se = I.se,
              rg = rg,
              rg.se = rg.se)
         )
}


