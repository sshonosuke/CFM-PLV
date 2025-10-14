###---------------------------------------------------------------------###
###       Code for Experiment: PLV estimation with noisy data           ###
###---------------------------------------------------------------------###
rm(list=ls())
library(MCMCpack)
library(progress)

## load R functions
source("WMF.R")
set.seed(1)

## settings 
error <- "Gaussian"   # "Gaussian" or "uniform"
noise_level <- 0.1    # 0.1, 0.2, 0.3, 0.4, 0.5, 0.6 for Gaussian / 0.2, 0.4, 0.6, 0.8, 1.0, 1.2 for Uniform
data_type <- "alpha"  # alpha or beta


## data loading
load(paste0("Data-", data_type, "-", error, "-noise", noise_level, ".RData"))
n <- dim(Y_n)[[1]]
p <- dim(Y_n)[[2]]
TT <- dim(Y_n)[[3]]

## analysis 
draw <- 1500    # number of posterior draws
burn <- 500     # number of burn-in draws
fit <- WMF(Y=Y_n, L=12, draw=draw, burn=burn, spline=T)


## PLV
mc <- draw - burn
sub <- upper.tri(matrix(NA, p, p))
PLV_ind_pos_mat <- array(NA, c(mc, n, p, p))
for(k in 1:p){
  for(j in 1:p){
    for(i in 1:n){
      for(r in 1:mc){
        th1 <- fit$theta[r,i,k,] - 2*pi*fit$Z[r,i,k,]
        th2 <- fit$theta[r,i,j,] - 2*pi*fit$Z[r,i,j,]
        PLV_ind_pos_mat[r,i,k,j] <- abs( mean(exp(1i*(th1-th2))) )
      }
    }
  }
}

PLV_pos <- apply(PLV_ind_pos_mat, c(1,3,4), mean)
hPLV_model <- c(apply(PLV_pos, c(2,3), mean)[sub])


# PLV (direct calculation without denoising) 
hPLV_naive <- matrix(NA, p, p)
for(k in 1:p){
  for(j in 1:p){
    hPLV_naive[k,j] <- mean( abs(apply(exp(1i*(Y_n[,k,]-Y_n[,j,])), 1, mean)) )
  }
}
hPLV_naive <- hPLV_naive[sub]


## Save 
save(PLV_pos, hPLV_model, hPLV_naive, file=paste0("Result-", data_type, "-", error, "-noise", noise_level, ".RData"))


