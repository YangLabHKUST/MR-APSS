# The function for calculating likelihood
cal_likeli <- function(hat.b, A, Var, S, C1, C2, pi0){

  S2 <- mapply(function(M1,M2){M1 +  A %*% M2 %*% t(A) }, M1=S, M2=Var)
  S2 = lapply(split(S2, col(S2)), function(x){matrix(x, nrow=2, ncol=2)})

  P1 =  mapply(mvtnorm::dmvnorm, x= hat.b,  sigma=S2)
  P2 =  mapply(mvtnorm::dmvnorm, x= hat.b,  sigma=S)

  #return likelihood
  return(sum(log(pi0 * P1/2/pnorm(C1) + (1-pi0) * P2/2/pnorm(C2))))

}


MRAPSS_EM_fixr_func <- function(data = NULL,
                                beta =NULL,
                                Sigma = NULL,
                                pi0 = NULL,
                                fix.beta=T,
                                fix.r = T,
                                r = 0,
                                C = diag(2),
                                Omega = matrix(0,2,2),
                                tol=1e-06,
                                Threshold=1,
                                ELBO=F){

  if(is.null(data)){
    message("No data for MR testing")
    return(NULL)
  }

  m=nrow(data)
  Cr <- abs(qnorm(Threshold/2))

  se.exp = data$se.exp
  se.out = data$se.out

  LDsc = data$L2
  if(is.null(LDsc)) LDsc = rep(1, m)
  LDsc = ifelse(LDsc<1, 1, LDsc)

  # genome-wide shared
  Omega = matrix(LDsc, nrow=m, ncol=1) %*% matrix(as.vector(Omega), nrow = 1)

  # error term covaraince matrix: s_j(\rho)
  s = matrix(0, m, 4)
  s[1:m, 1] = se.exp^2 * drop(C[1,1])
  s[1:m, 4] = se.out^2 * drop(C[2,2])
  s[1:m, 2] = s[1:m, 3] =  drop(C[1,2]) * se.exp * se.out
  s11 =  s[1:m, 1]

  # genome wide shared + s_j(\rho)
  S = matrix(s + Omega, nrow=m, ncol=4)
  S11 = S[,1]



  R = matrix(c(1, r, r, 1), nrow = 2, ncol = 2)
  H = t(chol(R))
  inv.H = solve(H)

  # initialize
  if(is.null(beta)){
    beta = 0
  }

  if(is.null(Sigma)){
    sigma.sq = drop(max(mean(data$b.exp^2) - mean(S[,1]), 1e-8))/mean(LDsc)
    tau.sq =  drop(max(mean(data$b.out^2) - mean(S[,4]), 1e-8))/mean(LDsc)
    Sigma = matrix(c(sigma.sq, r*sqrt(sigma.sq*tau.sq), r*sqrt(sigma.sq*tau.sq), tau.sq), 2,2)
  }

  if(is.null(pi0))   pi0 = 0.5

  V1 = matrix(c(0,0,1,0),2,2,byrow = T)

  S = lapply(split(S, row(S)), function(x){matrix(x, nrow=2, ncol=2)})

  inv.S <- lapply(S, solve)

  hat.b = t(as.matrix(data[,c("b.exp","b.out")]))

  hat.b = lapply(split(hat.b, col(hat.b)), function(x){c(x)})

  likelis = NULL
  elbos = NULL

  C_2 = - Cr*sqrt(s11)/sqrt(S11)

  for(i in 1:5000){
    sigma.sq = drop(Sigma[1,1])
    Var =  LDsc*matrix(as.vector(Sigma), nrow=m,ncol=4, byrow=T)

    inv.Var =  1/LDsc*matrix(as.vector(solve(Sigma)), nrow=m,ncol=4, byrow=T)

    Var = lapply(split(Var, row(Var)), function(x){matrix(x, nrow=2, ncol=2)})

    inv.Var = lapply(split(inv.Var, row(inv.Var)), function(x){matrix(x, nrow=2, ncol=2)})

    A = matrix(c(1, 0, beta, 1), 2,2, byrow = T)

    # E-step
    C_1 = - Cr*sqrt(s11)/sqrt(S11 + LDsc*sigma.sq)

    # E_p[(gamma_j, \alpha_j) | Z_j =1], Var_p[(gamma_j, \alpha_j) | Z_j =0]

    inv.Sigmaj <-  mapply(function(M1, M2){t(A) %*% M1 %*% A + M2}, M1=inv.S, M2=inv.Var)

    inv.Sigmaj <- lapply(split(inv.Sigmaj, col(inv.Sigmaj)), function(x){matrix(x, nrow=2, ncol=2)})

    Sigmaj <- lapply(inv.Sigmaj, solve)

    muj <- mapply(function(M1, M2, vec1){M2 %*% t(A) %*% M1 %*% vec1},
                  vec1=hat.b,
                  M1=inv.S,
                  M2=Sigmaj)

    muj = lapply(split(muj, col(muj)),
                 function(x){matrix(x, nrow=2, ncol=1)})

    # E_p(Z_j)
    if(pi0 == 1) Pi = rep(1, m)

    if(pi0 !=1)  {

      bj = 1/2 *mapply(function(M1, vec1){ t(vec1) %*% solve(M1) %*% vec1}, vec1 = muj, M1=Sigmaj) +
        log(pi0/(1-pi0)) + 1/2* mapply(function(M1, M2){ log(det(M1)) - log(det(M2))}, M1 = Sigmaj, M2 = Var) +
        log(pnorm(C_2)/pnorm(C_1))

      Pi = 1/(1+exp(-bj))

      Pi = ifelse(Pi<1e-04, 1e-04, Pi)
      Pi = ifelse(Pi>0.9999, 0.9999, Pi)

    }

    likeli <- cal_likeli(hat.b, A, Var, S, C_1, C_2, pi0)
    likelis <-  c(likelis, likeli)

    if(i>1 && likelis[i] < likelis[(i-1)]) {
      warning("Likelihood decreasing")
      cat("Likelihood: ", likeli, "\t","beta", beta, "\t", "pi", pi0, "\t" ,"r", r, "\n")
      print(Sigma)
      print("Likelihood decreasing")
    }
    if(i>1 && abs((likelis[i]-likelis[(i-1)])/likelis[(i-1)]) < tol)  break

    #cat("Marginal Likelihood: ", likeli, "\t", "ELBO:", elbo,"\t", "Difference: ",likeli-elbo ,"\n")
    #cat("beta", beta, "\t", "pi", pi0, "\t" ,"r", r, "\n")

    # M-step
    if(!fix.beta){
      # update beta
      temp1 = sum(Pi * mapply(function(vec1, M1,vec2)
        return(t(vec1)%*% t(V1) %*% M1 %*% vec2),
        vec1=muj, M1=inv.S, vec2= hat.b)) -
        sum(Pi * mapply(function(M1,M2,vec1)
          return(sum(diag(t(V1) %*% M1 %*% (M2 + vec1 %*% t(vec1))))),
          M1 =inv.S , M2 = Sigmaj, vec1 = muj))

      temp2 =
        sum(Pi * mapply(function(M1,M2,vec1)
          return(sum(diag(t(V1) %*% M1 %*% V1 %*% (M2 + vec1 %*% t(vec1))))),
          M1 =inv.S , M2 = Sigmaj, vec1 = muj))

      beta = temp1/temp2

    }


    ## update pi0 when pi0!=1
    if(pi0!=1){

      pi0 = sum(Pi)/m
      # set lower bound and upper bound of $\pi_0$ to avoid log(0)
      if(pi0<1e-04)   pi0 = 1e-04
      if(pi0>0.9999)  pi0 = 0.9999

    }

    ## update Sigma

    temp0 = mapply(function(vec1, M1) vec1 %*% t(vec1) + M1,
                   vec1 = muj, M1 = Sigmaj)

    tempL = matrix(temp0 %*% (Pi/LDsc), byrow=F, nrow = 2, ncol =2)

    if(Threshold==1) {

      Sigma = tempL/sum(Pi)
    }

    if(Threshold!=1) {

      E = matrix(0, 2, 2)
      E[1,1] = 1

      tempR = sum(Pi) * solve(Sigma) +  sum(Pi * dnorm(C_1)/pnorm(C_1) * Cr * LDsc *
                                              sqrt(s11) * (S11 + LDsc*sigma.sq)^(-3/2)) * E
      L = t(chol(tempR))
      inv.L = solve(L)

      Sigma = t(inv.L) %*% expm::sqrtm(t(L) %*% tempL %*% L) %*% inv.L

      if(fix.r & r==0){
        sigma.sq = sqrt(tempL[1,1]/tempR[1,1])
        tau.sq = tempL[2,2]/sum(Pi)
        Sigma = diag(c(sigma.sq, tau.sq))
      }

      if(fix.r & r!=0){
        #Sigma[1,2] = Sigma[2,1] = r*sqrt(Sigma[1,1]*Sigma[2,2])
        sqrt.Sigma = t(inv.H) %*% expm::sqrtm(t(H) %*% Sigma %*% H) %*% inv.H
        Sigma = diag(diag(sqrt.Sigma)) %*% R %*% diag(diag(sqrt.Sigma))
      }

  }
}



  return(list(beta = beta,
              Sigma = Sigma,
              pi0 = pi0,
              likeli=likeli,
              Threshold = Threshold))

}

