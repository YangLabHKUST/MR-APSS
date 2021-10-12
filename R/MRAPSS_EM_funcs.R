# The function for calculating evidence lower bound
cal_elbo <- function(A, pi0, m, Pi, S, inv.S, Var, inv.Var, muj, Sigmaj, hat.b, C1, C2){

  return(
    # Elog(p(\hat gamma_j, \hat Gamma_j | \gamma_j, \alpha_j, Z_j))
    -log(2*pi) * m - 1/2*sum(log(unlist(lapply(S, det)))) -
      1/2 *sum(mapply(function(M1, vec1){t(vec1) %*% M1 %*% vec1}, vec1=hat.b, M1=inv.S)) +
      sum(Pi *mapply(function(M1, vec1, vec2){t(vec1) %*% M1 %*% A %*% vec2}, vec1=hat.b, M1=inv.S, vec2=muj)) -
      1/2 * sum(Pi*mapply(function(M1, vec1, M2){sum(diag(t(A) %*% M1 %*% A %*% (vec1 %*% t(vec1) + M2)))},
                          vec1=muj, M1=inv.S, M2=Sigmaj)) +
      # E(log p((Z_j))
      sum(Pi * log(pi0 + (pi0==0)) + (1-Pi)*log(1-pi0 + (pi0==1))) +
      # E(log p(gamma_j, alpha_j))
      -log(2*pi) * m - 1/2*sum(log(unlist(lapply(Var, det)))) -
      1/2 * sum(Pi * mapply(function(M1, vec1, M2){sum(diag(M1 %*% (vec1 %*% t(vec1) + M2)))},
                            vec1=muj, M1=inv.Var, M2=Sigmaj)) -
      sum(1-Pi) +
      # -E(log(q(Z_j)))
      sum(-Pi * log(Pi + (Pi==0)) - (1-Pi) * log(1-Pi + (Pi==1))) +
      # -E(log(q(gamma_j,alpha_j|Z_j)))
      sum(Pi * (1 + log(2*pi) + 1/2 *log(unlist(lapply(Sigmaj, det))))) +
      sum((1-Pi) * (1 + log(2*pi) + 1/2 *log(unlist(lapply(Var, det))))) -
      sum(Pi *log(2*pnorm(C1)) + (1-Pi)*log(2*pnorm(C2)))
  )


}

MRAPSS_EM_func <- function(data = NULL,
                          beta=NULL,
                          sigma.sq = NULL,
                          tau.sq = NULL,
                          pi0 = NULL,
                          fix.beta = F,
                          fix.tau=F,
                          fix.sigma=F,
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


  # initialize
  if(is.null(beta)){
    beta = 0
  }

  if(is.null(sigma.sq) ){
    sigma.sq = drop(max(mean(data$b.exp^2) - mean(S[,1]), 1e-8))/mean(LDsc)
  }

  if(is.null(tau.sq)){
    tau.sq =  drop(max(mean(data$b.out^2) - mean(S[,4]), 1e-8))/mean(LDsc)
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

    Var =  LDsc*matrix(as.vector(diag(c(sigma.sq,tau.sq))), nrow=m,ncol=4, byrow=T)

    inv.Var =  1/LDsc*matrix(as.vector(diag(c(1/sigma.sq,1/tau.sq))), nrow=m,ncol=4, byrow=T)

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
    #if(!ELBO){
      if(i>1 && likelis[i] < likelis[(i-1)])  warning("Likelihood decreasing")
      if(i>1 && abs((likelis[i]-likelis[(i-1)])/likelis[(i-1)]) < tol)  break
    #}else{
      elbo <- cal_elbo(A, pi0, m, Pi, S, inv.S, Var, inv.Var, muj, Sigmaj, hat.b, C_1, C_2)
      elbos <-  c(elbos, elbo)
      #if(i>1 && elbos[i] < elbos[(i-1)])  warning("ELBO decreasing")
      #if(i>1 && abs((elbos[i]-elbos[(i-1)])/elbos[(i-1)]) < tol)  break
    #}

    #cat("Marginal Likelihood: ", likeli, "\t", "ELBO:", elbo,"\t", "Difference: ",likeli-elbo ,"\n")
    #cat("beta", beta, "\t", "pi", pi0, "\t" ,"sigma.sq", sigma.sq, "\t", "tau.sq", tau.sq, "\n")

    # M-step
    if(fix.beta==F){

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

    ## update sigma^2
    mu.tr = data.table::transpose(muj)
    mu.vec = as.matrix(data.frame(mu.tr[[1]], mu.tr[[2]]))

    Sigma.tr = data.table::transpose(Sigmaj)
    Sigma.vec = as.matrix(data.frame(Sigma.tr[[1]], Sigma.tr[[4]]))

    if(sigma.sq==0) fix.sigma=T
    if(tau.sq==0) fix.tau = T

    if(!fix.sigma){

      A1 = sum(Pi * mu.vec[,1]^2/LDsc + Pi*Sigma.vec[,1]/LDsc)

      if(Threshold==1)  sigma.sq = A1/sum(Pi)

      if(Threshold!=1) {

        D1 = sum(Pi)/sigma.sq

        B1 =  sum(Pi * dnorm(C_1)/pnorm(C_1) * Cr * LDsc *
                    sqrt(s11) * (S11 + LDsc*sigma.sq)^(-3/2))

        sigma.sq = sqrt(A1/(D1 + B1))

      }

    }

    # update tau^2
    if(!fix.tau){

      A2 = sum(Pi*mu.vec[,2]^2/LDsc + Pi*Sigma.vec[,2]/LDsc)

      tau.sq = A2/sum(Pi)

    }

  }

  return(list(beta = beta,
              sigma.sq=sigma.sq,
              tau.sq =tau.sq,
              pi0 = pi0,
              fix.tau = fix.tau,
              fix.sigma = fix.sigma,
              post = list(mu = mu.vec, Pi = Pi, IVsignal.sum =  sum(Pi * mu.vec[,1]^2 + Pi*Sigma.vec[,1])),
              likeli = likeli,
              likelis = likelis,
              elbos = elbos,
              Threshold = Threshold))

}

