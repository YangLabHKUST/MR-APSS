get.h2.weights <- function(h2, intercept, L2, N, M, x){
  #
  h2 <- max(h2, 0)
  h2 <- min(h2, 1)
  wld <- as.numeric(lapply(X=L2, function(x){max(x,1)}))
  ld <- as.numeric(lapply(X=x, function(x){max(x,1)}))
  #cat("ld:",mean(ld),"\n")
  c <- h2*N/M
  het.w <- 1/(2*(intercept + c*ld )^2)
  #cat("mean het.w ", mean(het.w),"\n")
  oc.w <- 1/wld
  #cat("mean oc.w ", mean(oc.w), "\n")
  w <- sqrt(het.w*oc.w)
  return(w)
}


irwls_h2 <- function(y, L2, update.x, weights, intercept=1,
                     M, N, N.bar, fix.intercept=T, seperator, jknife=T){

  # calculate new weights
  weights = weights/sum(weights)
  x = L2 *N/N.bar

  if(fix.intercept){
    y = (y-intercept)
    for(i in 1:2){

      wy = y*weights

      wx = as.matrix(x*weights)

      fit = lm(wy~wx+0)


      h2 = drop(coef(fit)[1]) * M/N.bar

      h2 <- max(h2, 0)

      h2 <- min(h2, 1)

      #cat(h2,"\t")
      #cat(intercept,"\n")
      # update weights
      weights = get.h2.weights(h2, intercept, L2, N, M, update.x)
      #cat(mean(weights), "\n")
      weights = weights/sum(weights)
      #print(intercept)
      #cat("i:", i," ", mean(weights), "\n")
    }

    ## preweight LD and chi:
    weighted.LD <- as.matrix(x*weights)
    weighted.chi <- as.matrix(y*weights)

  }else{

    for(i in 1:2){

      wy = y*weights

      wx = as.matrix(cbind(x*weights, weights))

      fit = lm(wy~wx+0)
      h2 = coef(fit)[1]*M/N.bar

      h2 <- max(h2, 0)
      h2 <- min(h2, 1)
      intercept = coef(fit)[2]
      #cat(h2,"\t")
      #cat(intercept,"\n")
      # update weights
      weights = get.h2.weights(h2, intercept, L2, N, M, update.x)
      #cat(mean(weights), "\n")
      weights = weights/sum(weights)
      #cat("i:", i," ", mean(weights), "\n")
    }


    ## preweight LD and chi:

    weighted.LD <- as.matrix(cbind(x*weights, weights))
    weighted.chi <- as.matrix(y*weights)

  }


  ## Perfrom analysis:

  n.blocks = length(seperator)-1
  n.snps <- length(y)
  p = ncol(weighted.LD)
  xty.block.values <- matrix(data=NA, nrow=n.blocks, ncol=p)
  xtx.block.values <- matrix(data=NA, nrow=p*n.blocks, ncol=p)

  from <- seperator
  to <- c(seperator[2:n.blocks]-1, n.snps)

  rep.from <- seq(from=1,to=p*n.blocks , by=p)
  rep.to <- seq(from =p,to=p*n.blocks ,by =p)

  colnames(xty.block.values)<- colnames(xtx.block.values)<- colnames(weighted.LD)


  for(i in 1:n.blocks){
    xty.block.values[i,] <- t(t(weighted.LD[from[i]:to[i],]) %*% weighted.chi[from[i]:to[i],])
    xtx.block.values[rep.from[i]:rep.to[i],] <- as.matrix(t(weighted.LD[from[i]:to[i],]) %*%
                                                            weighted.LD[from[i]:to[i],])
  }

  xty <- as.matrix(colSums(xty.block.values))
  xtx <- matrix(data=NA,nrow =p, ncol =p)
  colnames(xtx)<- colnames(weighted.LD)

  for(i in 1:p){
    xtx[i,] <- t(colSums(as.matrix(xtx.block.values[seq(from=i, to=p*n.blocks, by = p),])))
  }


  reg <- solve(xtx)%*% xty

  if(jknife){


    # perform jacknife
    delete.from <- seq(from=1,to= p* n.blocks, by=p)
    delete.to <- seq(from=p,  to=p* n.blocks, by=p)
    delete.values <- matrix(data=NA, nrow=n.blocks, ncol = p)
    colnames(delete.values)<- colnames(weighted.LD)

    for(i in 1:n.blocks){
      xty.delete <- xty-xty.block.values[i,]
      xtx.delete <- xtx-xtx.block.values[delete.from[i]:delete.to[i],]
      delete.values[i,] <- solve(xtx.delete)%*% xty.delete
    }

    delete.values <- as.matrix(delete.values[,1:p])

    pseudo.values <- matrix(data=NA,nrow=n.blocks,ncol=p)
    colnames(pseudo.values)<- colnames(weighted.LD)
    for(i in 1:n.blocks){
      pseudo.values[i,] <- (n.blocks*reg)-((n.blocks-1)* delete.values[i,])
    }

    jackknife.cov <- cov(pseudo.values)/n.blocks
    jackknife.se <- sqrt(diag(jackknife.cov))

    coef.cov <- jackknife.cov[1,1]/(N.bar^2)*M^2
    h2.se <- sqrt(coef.cov)

    if(p==1){
      intercept.est <- intercept
      h2.est <- reg[1]/N.bar*M
      intercept.se <- NA
    }

    if(p==2){

      intercept.est <- reg[2]
      h2.est <- reg[1]/N.bar*M
      intercept.se <- jackknife.se[length(jackknife.se)]

    }

    return(list(h2 = h2.est, h2.se= h2.se,
                intercept = intercept.est,
                intercept.se = intercept.se,
                delete.values =  delete.values,
                jk.est = reg[1]))
  }else{


    if(p==1){
      intercept.est <- intercept
      h2.est <- reg[1]/N.bar*M
    }

    if(p==2){

      intercept.est <- reg[2]
      h2.est <- reg[1]/N.bar*M

    }
    return(list(h2 = h2.est,
                intercept = intercept.est))

  }

}


est_h2 <- function(merged, M, trait.name = NULL, Twostep=F, fix.intercept=T, step1.idx=NULL, n.blocks = 200, jknife=T){

  ## Read in summary statistics
  require(readr)
  if(is.null(trait.name)) trait.name = Trait
  colnames(merged) = c("SNP","chi2","N","L2")
  chi2 = merged$chi2
  N = merged$N
  L2 = merged$L2
  n.snps = nrow(merged)
  L2 <- as.numeric(lapply(X=L2, function(x){max(x,1)}))

  ## initial WEIGHTS:
  intercept = 1
  tot.agg <- (M*(mean(chi2)-1))/mean(L2*N)
  weights = get.h2.weights(tot.agg, intercept = 1, L2, N, M,L2)
  N.bar <- mean(N)
  x = L2 *N/N.bar

  if(fix.intercept) Twostep = F
  if(is.null(step1.idx)) Twostep = F
  if(Twostep){
    ### Two dtep estimator
    ## step I
    chi2.s1 <- chi2[step1.idx]
    L2.s1   <- L2[step1.idx]
    weights.s1 <- weights[step1.idx]
    N.s1 <- N[step1.idx]
    x.s1 <- x[step1.idx]
    seperator <- floor(seq(from=1, to=length(step1.idx), length.out =(n.blocks+1)))
    step1 = irwls_h2(chi2.s1, L2 = L2.s1, update.x=x.s1, weights.s1, intercept=1, M, N.s1, N.bar, fix.intercept, seperator, jknife)

    ## step 2
    new.seperator <-  c(1, step1.idx[seperator[2:n.blocks]], n.snps)
    step2 = irwls_h2(chi2,  L2=L2, update.x=L2, weights, intercept=step1$intercept, M, N, N.bar, fix.intercept=T, new.seperator,jknife)
    h2 = step2$h2
    intercept = step1$intercept

    if(jknife){
      ## cpmbine step 1 and step 2
      c = sum(weights^2 *x)/sum(weights^2*x^2)

      est = c(step2$jk.est,step1$intercept)
      delete.values = matrix(0, nrow = n.blocks, ncol = 2)
      delete.values[,2] = step1$delete.values[,2]
      delete.values[,1] = (step2$delete.values - c * (step1$delete.values[,2] - step1$intercept))

      pseudo.values = matrix(n.blocks*est,nrow =n.blocks,ncol=2, byrow = T) - (n.blocks-1)* delete.values


      jackknife.cov <- cov(pseudo.values)/n.blocks
      jackknife.se <- sqrt(diag(jackknife.cov))

      coef.cov <- jackknife.cov[1,1]/(N.bar^2)*M^2
      h2.se <- sqrt(coef.cov)
      intercept.se <- jackknife.se[2]

    }

  }else{
    seperator <- floor(seq(from=1, to=length(chi2), length.out =(n.blocks+1)))
    step1 = irwls_h2(chi2,  L2=L2, update.x=L2, weights, intercept=1, M, N, N.bar, fix.intercept, seperator)
    h2 = step1$h2
    intercept = step1$intercept
    if(jknife){
      h2.se = step1$h2.se
      intercept.se = step1$intercept.se
      delete.values = step1$delete.values
    }

  }


  # lambda.gc <- median(merged$chi2)/0.4549
  # mean.Chi <- mean(merged$chi2)
  # ratio <- (intercept-1)/(mean(merged$chi2)-1)
  # ratio.se <- intercept.se/(mean(merged$chi2)-1)

  if(jknife){
    return(list(h2 = h2, h2.se = h2.se,
                intercept = intercept, intercept.se =intercept.se,
                delete.values= delete.values[,1]))
  }else{
    return(list(h2 = h2,
                intercept = intercept))
  }



}
