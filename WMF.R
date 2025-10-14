###---------------------------------------------------------------------###
###        Code for Wrapped Multivariate Functional Model               ###
###---------------------------------------------------------------------###
library(MCMCpack)
library(progress)
library(splines)
library(signal)

WMF <- function(Y, L=10, draw=1500, burn=500, Basis=NULL, spline=F){
  ## preparation 
  n <- dim(Y)[1]
  p <- dim(Y)[2]
  TT <- dim(Y)[3]
  MC <- draw - burn 
  if(is.null(Basis)){
    if(spline){  # quadratic spline
      t_vec <- (1:TT)/TT
      B <- cbind(1, t_vec, t_vec^2)
      qq <- dim(B)[2]
      K <- L-qq    # number of knots
      knot <- seq(0, 1, length=K+2)[2:(K+1)]
      for(k in 1:K){
        B <- cbind(B, (t_vec-knot[k])^2*(t_vec>knot[k]))
      }
    }else{
      B <- bs(1:TT, df=L, degree=3, intercept=T)   # B-spline
    }
  }else{
    B <- Basis
  }
  L <- dim(B)[2]
  zk <- seq(-3, 3, by=1)   # candidate integer for Z
  M <- length(zk)
  B_mat <- t(B)%*%B
  
  # prior 
  nu0 <- eta0 <- 1
  A0 <- 0
  B0 <- 10^4
  
  # initial value
  mu <- matrix(0, p, L)
  tau <- rep(1, L)    # std of individual/dimensional effects 
  gam <- rep(1, L)    # std of dimensional effects
  beta <- rep(0, L)   # global mean of coefficient
  sig <- 0.2
  Z <- array(0, c(n, p, TT))
  lam <- 1    # only for spline 
  
  # initial value of A 
  A <- array(0, c(n, p, L))
  wY <- Y
  for(i in 1:n){
    for(k in 1:p){
      wY[i,k,] <- unwrap(Y[i,k,])
    }
  }
  if(TT>L){
    for(i in 1:n){
      for(k in 1:p){
        A[i,k,] <- coef(lm(wY[i,k,]~B-1))
      }
    }
  }
  
  # objects for posterior samples
  A_pos <- array(NA, c(MC, n, p, L))
  mu_pos <- array(NA, c(MC, p, L))
  tau_pos <- gam_pos <- beta_pos <- matrix(NA, MC, L)
  sig_pos <- rep(NA, MC)
  theta_pos <- array(NA, c(MC, n, p, TT))
  Z_pos <- array(NA, c(MC, n, p, TT))
  
  ## MCMC iteration
  pb <- progress_bar$new(total=draw)   # progress bar 
  for(itr in 1:draw){
    # mu
    sA <- apply(A, c(2,3), sum)
    for(l in 1:L){
      pos_var_l <- 1 / (1/gam[l]^2 + n/tau[l]^2)
      pos_mean_l <- pos_var_l * (beta[l]/gam[l]^2 + sA[,l]/tau[l]^2)
      mu[,l] <- rnorm(p, pos_mean_l, sqrt(pos_var_l))
    }
    
    # beta 
    prior_var <- rep(B0, L)
    if(spline){
      prior_var[-(1:qq)] <- lam
    }
    pos_var <- 1/ (p/gam^2 + 1/prior_var) 
    pos_mean <- pos_var * (apply(mu, 2, sum)/gam^2 + A0/prior_var)
    beta <- rnorm(L, pos_mean, sqrt(pos_var))
    if(spline){
      lam <- rinvgamma(1, 2+0.5*K, 0.1+0.5*sum(beta[-(1:qq)]^2))
    }
    
    # tau 
    nu1 <- nu0 + 0.5*n*p
    for(l in 1:L){
      eta1 <- eta0 + 0.5*sum((t(A[,,l])-mu[,l])^2)
      tau[l] <- sqrt( rinvgamma(1, nu1, eta1) )
    }
    
    # gamma 
    nu1 <- nu0 + 0.5*p
    eta1 <- eta0 + 0.5*apply((t(mu)-beta)^2, 1, sum)
    gam <- sqrt( rinvgamma(L, nu1, eta1) )
    
    # Z
    for(i in 1:n){
      rr <- Y[i,,] - A[i,,]%*%t(B)
      shift_int <- (-1)*round(rr/(2*pi))  # (p,TT)-matrix for the i-th obs
      dens <- array(NA, c(p, TT, M))
      for(m in 1:M){
        dens[,,m] <- dnorm(rr+2*pi*(shift_int+zk[m]), 0, sig)
      }
      for(k in 1:p){
        denom <- apply(dens[k,,], 1, sum)
        prob <- dens[k,,]/denom
        Z[i,k,] <- shift_int[k,] + apply(prob, 1, sample, x=zk, size=1, replace=F)
      }
    }
    
    # A
    At <- solve( B_mat/sig^2 + diag(1/tau^2) )
    for(k in 1:p){
      for(i in 1:n){
        Bt <- as.vector( t(B)%*%(Y[i,k,]+2*pi*Z[i,k,]) )/sig^2 + mu[k,]/tau^2
        A[i,k,] <- mvrnorm(1, At%*%Bt, At)
      }
    }
    
    # regression term 
    theta <- array(NA, c(n, p, TT))
    for(i in 1:n){
      theta[i,,] <- A[i,,]%*%t(B)
    }
    
    # sigma
    nu1 <- nu0 + n*p*TT/2
    eta1 <- eta0 
    for(i in 1:n){
      resid <- Y[i,,] - A[i,,]%*%t(B) + 2*pi*Z[i,,]
      eta1 <- eta1 + sum(resid^2)/2
    }
    sig <- sqrt( rinvgamma(1, nu1, eta1) )
    
    # save
    if(itr>burn){
      cc <- itr - burn
      mu_pos[cc,,] <- mu
      tau_pos[cc,] <- tau
      gam_pos[cc,] <- gam
      beta_pos[cc,] <- beta
      A_pos[cc,,,] <- A
      sig_pos[cc] <- sig
      theta_pos[cc,,,] <- theta
      Z_pos[cc,,,] <- Z
    }
    pb$tick()
  }
  
  # output 
  Result <- list(theta=theta_pos, A=A_pos, tau=tau_pos, mu=mu_pos, 
                 sig=sig_pos, Base=B, Z=Z_pos)
  return(Result)
}

