get.gencov.weights <- function(L2, h1, h2, intercept.h1,
                               intercept.h2, intercept.gencov,
                               N.x, N.y, rho_g, N.bar, M, x){
  #### MAKE WEIGHTS:
  h1 <- max(h1, 0)
  h1 <- min(h1, 1)

  h2 <- max(h2, 0)
  h2 <- min(h2, 1)

  rho_g <- max(rho_g, -1)
  rho_g <- min(rho_g, 1)
  wld <- as.numeric(lapply(X=L2, function(x){max(x,1)}))
  ld <- as.numeric(lapply(X=x,function(x){max(x,1)}))
  a = N.x * h1 * ld/M + intercept.h1
  b = N.y * h2 * ld/M + intercept.h2
  c = sqrt(N.x*N.y) * rho_g * ld/M + intercept.gencov
  het.w <- 1/(a*b + c^2)
  oc.w <- 1/wld
  #cat("mean het.w ", mean(het.w),"\n")
  #cat("mean oc.w ", mean(oc.w), "\n")
  w <- sqrt(het.w * oc.w )

  return(w)

}

irwls_gencov <- function(y, L2, update.x, N,
                         weights,
                         N.x,N.y,
                         h1, h2,
                         intercept.h1,
                         intercept.h2,
                         intercept.gencov=0,
                         M, N.bar,
                         fix.intercept=T,
                         seperator,
                         jknife=T){
  #cat("Weights \n")
  #cat(mean(weights),"\n")
  weights = weights/sum(weights)

  x = L2 *N/N.bar

  if(fix.intercept){

    y= y - intercept.gencov

    for(i in 1:2){

      wy = y*weights

      wx = as.matrix(x*weights)

      fit = lm(wy~wx+0)

      rho_g = drop(coef(fit)[1]* M/N.bar)

      rho_g <- max(rho_g, -1)
      rho_g <- min(rho_g, 1)

      # update weights
      weights = get.gencov.weights(L2, h1, h2,
                                   intercept.h1,
                                   intercept.h2,
                                   intercept.gencov,
                                   N.x, N.y,
                                   rho_g,
                                   N.bar, M, update.x)
      #cat(mean(weights),"\n")
      weights <-  weights/sum(weights)
    }

    ## preweight LD and chi:
    weighted.LD <- as.matrix(x*weights)

    weighted.chi <- as.matrix(y*weights)



  }else{


    for(i in 1:2){

      wy = y*weights

      wx = as.matrix(cbind(x*weights, 1*weights))

      fit = lm(wy~wx+0)

      rho_g = drop(coef(fit)[1]* M/N.bar)

      rho_g <- max(rho_g, -1)
      rho_g <- min(rho_g, 1)

      intercept.gencov = drop(coef(fit)[2])

      weights = get.gencov.weights(L2,h1, h2,
                                   intercept.h1,
                                   intercept.h2,
                                   intercept.gencov,
                                   N.x, N.y,
                                   rho_g,
                                   N.bar, M, update.x)
      #cat(mean(weights),"\n")
      weights <-  weights/sum(weights)
    }
    ## preweight LD and chi:

    weighted.LD <- as.matrix(cbind(x*weights, 1*weights))

    weighted.chi <- as.matrix(y*weights)
  }

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
    rho_g.se <- sqrt(coef.cov)

    if(p==1){
      intercept.est <- intercept.gencov
      rho_g.est <- reg[1]/N.bar*M
      intercept.se <- NA
    }

    if(p==2){

      intercept.est <- reg[2]
      rho_g.est <- reg[1]/N.bar*M
      intercept.se <- jackknife.se[length(jackknife.se)]

    }
    return(list(rho_g = rho_g.est, rho_g.se= rho_g.se,
                intercept = intercept.est,
                intercept.se = intercept.se,
                delete.values =  delete.values,
                jk.est = reg[1]))
  }else{
    if(p==1){
      intercept.est <- intercept.gencov
      rho_g.est <- reg[1]/N.bar*M
    }

    if(p==2){

      intercept.est <- reg[2]
      rho_g.est <- reg[1]/N.bar*M
    }
    return(list(rho_g = rho_g.est,
                intercept = intercept.est))
  }

}

est_gencov <- function(merged,
                       h1, h2,
                       intercept.h1,
                       intercept.h2,
                       intercept.gencov=0,
                       M,Twostep=F,
                       fix.gcov.intercept=F,
                       step1.idx=NULL,
                       n.blocks=200,
                       jknife=T){

  Zxy = merged$Zxy
  L2 = merged$L2
  N.x = merged$N.x
  N.y = merged$N.y
  N = sqrt(N.x * N.y)
  N.bar <- mean(N)
  x = L2 *N/N.bar
  rho_g <- mean(Zxy)*M/mean(N*L2)

  weights = get.gencov.weights(L2, h1, h2,
                               intercept.h1,
                               intercept.h2,
                               intercept.gencov,
                               N.x, N.y,
                               rho_g, N.bar, M, L2)

  if(fix.gcov.intercept)  Twostep=F
  if(is.null(step1.idx)) Twostep = F
  if(Twostep==T){
    message("Using two-step estimator with cutoff at 30.")
    ### Two step estimator
    Zxy.s1 <-Zxy[step1.idx]
    L2.s1   <- L2[step1.idx]
    weights.s1 <- weights[step1.idx]
    N.s1 <- N[step1.idx]
    x.s1 <- x[step1.idx]
    N.x.s1 <- N.x[step1.idx]
    N.y.s1 <- N.y[step1.idx]
    n.snps <- nrow(merged)
    seperator <- floor(seq(from=1, to=length(step1.idx), length.out =(n.blocks+1)))

    step1 = irwls_gencov( y = Zxy.s1,
                          L2 = L2.s1,
                          update.x = x.s1,
                          N=N.s1,
                          weights = weights.s1,
                          N.x = N.x.s1,
                          N.y = N.y.s1,
                          h1, h2,
                          intercept.h1,
                          intercept.h2,
                          intercept.gencov=0,
                          M, N.bar,
                          fix.intercept=F,
                          seperator)


    ## step 2
    #cat("step II \n")
    new.seperator <-  c(1, step1.idx[seperator[2:n.blocks]], n.snps)

    step2 =irwls_gencov(y=Zxy,
                        L2=L2,
                        update.x = L2,
                        N = N,
                        weights = weights,
                        N.x = N.x,
                        N.y = N.y,
                        h1, h2,
                        intercept.h1,
                        intercept.h2,
                        intercept.gencov=step1$intercept,
                        M, N.bar,
                        fix.intercept=T,
                        new.seperator)
    rho_g = step2$rho_g
    intercept = step1$intercept

    if(jknife){
      ## cpmbine step 1 and step 2
      x = L2 *N/N.bar
      c = sum(weights^2 *x)/sum(weights^2*x^2)
      est = c(step2$jk.est, step1$intercept)
      delete.values = matrix(0, nrow = n.blocks, ncol = 2)
      delete.values[,2] = step1$delete.values[,2]
      delete.values[,1] = M/N.bar *(step2$delete.values - c * (step1$delete.values[,2] - step1$intercept))

      pseudo.values = matrix(n.blocks*est,nrow =n.blocks,ncol=2, byrow = T) - (n.blocks-1)* delete.values

      jackknife.cov <- cov(pseudo.values)/n.blocks
      jackknife.se <- sqrt(diag(jackknife.cov))

      coef.cov <- jackknife.cov[1,1]
      rho_g.se <- sqrt(coef.cov)
      intercept.se <- jackknife.se[2]

    }

  }else{

    seperator <- floor(seq(from=1, to=nrow(merged), length.out =(n.blocks+1)))

    step1 = irwls_gencov( y = Zxy,
                          L2 = L2,
                          update.x = L2,
                          N=N,
                          weights = weights,
                          N.x = N.x,
                          N.y = N.y,
                          h1, h2,
                          intercept.h1,
                          intercept.h2,
                          intercept.gencov,
                          M, N.bar,
                          fix.intercept=fix.gcov.intercept,
                          seperator)
    rho_g = step1$rho_g
    intercept = step1$intercept
    if(jknife){
      rho_g.se = step1$rho_g.se
      intercept.se = step1$intercept.se
      delete.values=step1$delete.values
    }

  }

  if(jknife){
    return(list(rho_g = rho_g,
                rho_g.se = rho_g.se,
                intercept = intercept,
                intercept.se = intercept.se,
                delete.values=delete.values))
  }else{
    return(list(rho_g = rho_g,
                intercept = intercept))
  }


}


