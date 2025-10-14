###---------------------------------------------------------------------###
###           Code for summarizing experimental results                 ###
###---------------------------------------------------------------------###
rm(list=ls())
error <- "Gaussian"     # "Gaussian" or "uniform"
data_type <- "alpha"    #  "alpha" or "beta"


## baseline data (without noise)
load("PLV-original.RData")
if(data_type=="alpha"){ base <- hPLV_alpha }
if(data_type=="beta"){ base <- hPLV_beta }


## Absolute error 
M <- 6
AE <- CI95l <- CI95u <- matrix(NA, M, 2)

noise_level <- 0.1    # 0.1, 0.2, 0.3, 0.4, 0.5, 0.6 for Gaussian / 0.2, 0.4, 0.6, 0.8, 1.0, 1.2 for Uniform

if(erorr=="Gaussian"){ noise_level_set <- (1:M)/10 }
if(erorr=="uniform"){ noise_level_set <- 2*(1:M)/10 }
for(j in 1:M){
  load(paste0("Result-", data_type, "-", error, "-noise", noise_level_set[j], ".RData"))
  Est <- cbind(hPLV_model, hPLV_naive)
  AE[j,] <- apply(abs(Est-base), 2, mean)
  CI95l[j,] <- apply(abs(Est-base), 2, quantile, prob=0.05)
  CI95u[j,] <- apply(abs(Est-base), 2, quantile, prob=0.95)
}



## Plot
library(ggplot2)
library(dplyr)
library(tidyr)

method <- c("Model-based", "Naive")
df_plot <- data.frame(sigma = rep(sigma_set, times=2), method = rep(method, each=length(sigma_set)),
  AE = as.vector(AE), CI_low = as.vector(CI95l), CI_up = as.vector(CI95u))

# figure
pdf(paste0("Fig-", data_type, "-Gaussian.pdf"), height=7, width=9)
ggplot(df_plot, aes(x=sigma, y=AE, color=method, shape=method)) +
  geom_line(size=1) +
  geom_point(size=3) +
  geom_ribbon(aes(ymin=CI_low, ymax=CI_up, fill=method), alpha=0.2, color=NA) +
  labs(x="(Gaussian noise) band", y="Absolute Difference of PLV", color="Method", shape="Method", fill="Method") +
  theme_minimal(base_size=18) +
  scale_color_manual(values=c("blue", "red")) +
  scale_fill_manual(values=c("blue", "red")) +
  theme(legend.position="top")
dev.off()



## Threshold and posterior probability
K <- 0.7   # threshold
pos_prob_bound <- 0.5   # lower bound of posterior probability
p <- 90
sub <- upper.tri(matrix(NA, p, p))
Ind_base <- ifelse(base>K, 1, 0)   # base values

ER <- PLV <- TPR <- F1 <- matrix(NA, M, 2)
post_prob <- matrix(NA, sum(sub), M)
for(j in 1:M){
  load(paste0("Result-", data_type, "-", error, "-noise", noise_level_set[j], ".RData"))
  post_prob[,j] <- apply(PLV_pos>K, c(2,3), mean)[sub]
  Ind_pp <- ifelse(post_prob[,j]>pos_prob_bound, 1, 0)
  Ind <- cbind(Ind_pp, ifelse(hPLV_naive>K, 1, 0))
  ER[j,] <- 100*apply(abs(Ind-Ind_base), 2, mean)
  PLV[j,] <- 100*apply(Ind, 2, mean)
  TPR[j,] <- apply(Ind==1 & cbind(Ind_base, Ind_base)==1, 2, sum) / sum(Ind_base==1)
  precision <- apply(Ind==1 & cbind(Ind_base, Ind_base)==1, 2, sum) / apply(Ind==1, 2, sum)
  F1[j,] <- 2 * precision * TPR[j,] / (precision + TPR[j,])
}
F1[is.nan(F1)] <- 0

# table
write.csv(100*cbind(TPR, F1), file=paste0("Threshold-", data_type, ".csv"))



## reliability diagram (for alpha)
n_bins <- 5
bin_breaks <- seq(0, 1, length.out = n_bins + 1)

post_df <- as.data.frame(post_prob)
names(post_df) <- paste0(sigma_set)
post_df$y <- Ind_base
post_long <- post_df %>% pivot_longer(cols=-y, names_to="model", values_to="prob")

reliability_data <- post_long %>%
  mutate(bin = cut(prob, breaks = bin_breaks, include.lowest = TRUE)) %>%
  group_by(model, bin) %>%
  summarise(bin_mid = mean(prob), accuracy = mean(y), n=n(), .groups = "drop")

# figure
pdf(paste0("Reliability-", error, ".pdf"), height=7, width=7, pointsize=16)
ggplot(reliability_data, aes(x=bin_mid, y=accuracy, color=model)) +
  geom_line(size=1) +
  geom_point(size=2) +  
  geom_abline(intercept=0, slope=1, linetype="dashed", color="gray") +
  scale_x_continuous(limits=c(0,1), breaks=seq(0,1,1/n_bins)) +
  scale_y_continuous(limits=c(0,1), breaks=seq(0,1,1/n_bins)) +
  labs(
    x="Mean predicted probability",
    y="Observed accuracy",
    title="Gaussian noise",
    color="noise level",
    size="Bin count"
  ) +
  theme_minimal(base_size=14) +
  theme(plot.title=element_text(hjust=0.5))
dev.off()
