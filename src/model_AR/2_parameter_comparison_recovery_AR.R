## ---- PARAMETER COMPARISON RECOVERY: MODEL AR (ALPHA & RHO) ----
install.packages("pacman")
set.seed(42)
pacman::p_load(R2jags, parallel, polspline, tidyverse,ggplot2, 
               extraDistr, ggpubr, truncnorm, cascsim)
setwd("/work/LauraGivskovRahbek#6659/decision_making_exam24/src/model_AR")
source("../sim_functions.R")
source("../utils.R")

## ---- Simulation Specifications ----
ntokens <- 10 # represents the 10,000
vals <- seq(0, ntokens, 1) # values possible to contribute
groupsize <- 6 # number of participants per group
ngroups <- 4 # for one treatment (e.g. sunni in same-sect groups)
nsubs <- as.integer(groupsize*ngroups)
nsubs_m <- nsubs #two groups of equal size 
nsubs_s <- nsubs
niter <- 100

## ---- Empty Arrays for Simulated Data and True and Estimated Parameters ----
true_mu_alpha <- array(NA, c(niter))
true_delta_alpha <- array(NA, c(niter))
true_mu_rho <- array(NA, c(niter))
true_delta_rho <- array(NA, c(niter))
true_sigma_alpha_s <- array(NA, c(niter))
true_sigma_alpha_m <- array(NA, c(niter))
true_sigma_rho_s <- array(NA, c(niter))
true_sigma_rho_m <- array(NA, c(niter))

est_mu_alpha <- array(NA, c(niter))
est_delta_alpha <- array(NA, c(niter))
est_mu_rho <- array(NA, c(niter))
est_delta_rho <- array(NA, c(niter))
est_sigma_alpha_s <- array(NA, c(niter))
est_sigma_alpha_m <- array(NA, c(niter))
est_sigma_rho_s <- array(NA, c(niter))
est_sigma_rho_m <- array(NA, c(niter))

c_sim <- array(NA, c(2, niter, nsubs))
P_sim <- array(NA, c(2, niter, nsubs))
Gb_sim <- array(NA, c(2, niter, nsubs))
alpha_sim <- array(NA, c(2, niter, nsubs))
rho_sim <- array(NA, c(2, niter, nsubs))

## ---- Run n iterations of Recovery ----
for (i in 1:niter){
  
  # Sample mu (and delta) and sigma for ALPHA
  mu_alpha <- runif(1, min = 0.001, max = 10)
  delta_alpha <- rnorm(1, mean = 0, sd = 4)
  
  sigma_alpha_s <- rtruncnorm(1, a=0.001, b=Inf, mean = 1, sd = 2)
  sigma_alpha_m <- rtruncnorm(1, a=0.001, b=Inf, mean = 1, sd = 2)
  
  # Sample mu (and delta) and sigma for RHO
  mu_rho <- runif(1, min = 0.001, max = 0.999) #avoid extreme values (0 and 1)
  delta_rho <- rnorm(1, mean = 0, sd = 0.5) 
  sigma_rho_s <- runif(1, min = 10, max = 100)
  sigma_rho_m <- runif(1, min = 10, max = 100)
  
  # Simulate data based on parameters for group:s
  mu_alpha_s <- mu_alpha + (delta_alpha/2)
  mu_alpha_s <- ifelse(mu_alpha_s <=0.001, 0.001, ifelse(mu_alpha_s>10, 10, mu_alpha_s)) #avoid neg values by =0 if <0
  
  mu_rho_s <- mu_rho + (delta_rho/2)
  mu_rho_s <- ifelse(mu_rho_s <0, 0, ifelse(mu_rho_s>1, 1, mu_rho_s))
  
  sim_s <- PGGsimAR(nsubs_s, mu_alpha_s, mu_rho_s, sigma_alpha_s, sigma_rho_s)
  c_s <- sim_s$c
  
  # Save all simulations vars for inspection
  c_sim[1,i, ] <- c_s
  P_sim[1,i, ] <- sim_s$P
  Gb_sim[1,i, ] <- sim_s$Gb
  alpha_sim[1,i, ] <- sim_s$alpha
  rho_sim[1,i, ] <- sim_s$rho
  
  # Simulate data based on parameters for group:m
  mu_alpha_m <- mu_alpha - (delta_alpha/2)
  mu_alpha_m <- ifelse(mu_alpha_m <=0.001, 0.001, ifelse(mu_alpha_m>10, 10, mu_alpha_m))
  
  mu_rho_m <- mu_rho - (delta_rho/2)
  mu_rho_m <- ifelse(mu_rho_m <0, 0, ifelse(mu_rho_m >1, 1, mu_rho_m))
  
  sim_m <- PGGsimAR(nsubs_m, mu_alpha_m, mu_rho_m, sigma_alpha_m, sigma_rho_m)
  c_m <- sim_m$c
  
  # Save all simulations vars for inspection
  c_sim[2,i, ] <- c_m
  P_sim[2,i, ] <- sim_m$P
  Gb_sim[2,i, ] <- sim_m$Gb
  alpha_sim[2,i, ] <- sim_m$alpha
  rho_sim[2,i, ] <- sim_m$rho
  
  # Prep and run JAGS model
  data <- list("nsubs_s", "nsubs_m", "c_s", "c_m") 
  params <- c("mu_alpha", "mu_rho", "delta_alpha", "delta_rho",
              "sigma_alpha_s", "sigma_alpha_m", "sigma_rho_s", "sigma_rho_m")
  
  samples <- jags.parallel(data = data, inits=NULL, parameters.to.save = params, 
                           model.file ="PGGcompare_AR.txt", n.chains=3, n.iter=5000, 
                           n.burnin=1000, n.thin=1, jags.seed=3)
  
  # Save true and estimated parameters
  true_mu_alpha[i] <- mu_alpha
  true_delta_alpha[i] <- delta_alpha
  true_mu_rho[i] <- mu_rho
  true_delta_rho[i] <-delta_rho
  true_sigma_alpha_s[i] <- sigma_alpha_s
  true_sigma_alpha_m[i] <- sigma_alpha_m
  true_sigma_rho_s[i] <- sigma_rho_s
  true_sigma_rho_m[i] <- sigma_rho_m
  
  est_mu_alpha[i] <- MPD(samples$BUGSoutput$sims.list$mu_alpha)
  est_delta_alpha[i] <- MPD(samples$BUGSoutput$sims.list$delta_alpha)
  est_mu_rho[i] <- MPD(samples$BUGSoutput$sims.list$mu_rho)
  est_delta_rho[i] <- MPD(samples$BUGSoutput$sims.list$delta_rho)
  est_sigma_alpha_s[i] <- MPD(samples$BUGSoutput$sims.list$sigma_alpha_s)
  est_sigma_alpha_m[i] <- MPD(samples$BUGSoutput$sims.list$sigma_alpha_m)
  est_sigma_rho_s[i] <- MPD(samples$BUGSoutput$sims.list$sigma_rho_s)
  est_sigma_rho_m[i] <- MPD(samples$BUGSoutput$sims.list$sigma_rho_m)
  print(i)
}

## ---- Define vars for saving Data and Plotting Recovery ----
rdata_file <- "../../out/RData/AR_recov_data_comp.RData"
plot_path <- "../../out/AR_recovery_plot_comp.pdf"

param_names <- c("mu_alpha", "delta_alpha", "mu_rho", "delta_rho",
                 "sigma_alpha_s", "sigma_alpha_m", "sigma_rho_s", "sigma_rho_m") #names of parameters modeled

recov_titles <- c("Recovery of alpha mode (mu)", "Recovery of alpha treatment difference (delta)", 
                  "Recovery of rho mean (mu)", "Recovery of rho treatment difference (delta)",
                  "Recovery of alpha standard deviation (sigma) - same", 
                  "Recovery of alpha standard deviation (sigma) - mixed", 
                  "Recovery of rho concentration (sigma) - same",
                  "Recovery of rho concentration (sigma) - mixed")

## ---- Save Data ----
save(true_mu_alpha, true_delta_alpha, true_mu_rho, true_delta_rho, 
     true_sigma_alpha_s, true_sigma_alpha_m, true_sigma_rho_s, true_sigma_rho_m, 
     est_mu_alpha, est_delta_alpha, est_mu_rho, est_delta_rho, 
     est_sigma_alpha_s, est_sigma_alpha_m, est_sigma_rho_s, est_sigma_rho_m, 
     c_sim, P_sim, Gb_sim, alpha_sim, rho_sim,
     file = rdata_file)

## ---- Create Plots ----

names <- load(rdata_file, verbose = T) #names of parameter arrays 
param_df <- data.frame("true_mu_alpha" = true_mu_alpha, "est_mu_alpha" = est_mu_alpha, 
                       "true_delta_alpha" = true_delta_alpha, "est_delta_alpha" = est_delta_alpha,
                        "true_mu_rho" = true_mu_rho, "est_mu_rho" = est_mu_rho, 
                       "true_delta_rho" = true_delta_rho, "est_delta_rho" = est_delta_rho,
                        "true_sigma_alpha_s" = true_sigma_alpha_s, "est_sigma_alpha_s" = est_sigma_alpha_s,
                       "true_sigma_alpha_m" = true_sigma_alpha_m, "est_sigma_alpha_m" = est_sigma_alpha_m,
                       "true_sigma_rho_s" = true_sigma_rho_s, "est_sigma_rho_s" = est_sigma_rho_s,
                       "true_sigma_rho_m" = true_sigma_rho_m, "est_sigma_rho_m" = est_sigma_rho_m)

recov_plots <- list() #list to save recovery plots into
#min and max values for plots
min <- c(0, #mu_alpha
         -12, #delta_alpha
         0, #mu_rho
         -2, #delta_rho
         min(true_sigma_alpha_s, est_sigma_alpha_s), min(true_sigma_alpha_m, est_sigma_alpha_m), 
         min(true_sigma_rho_s, est_sigma_rho_s), min(true_sigma_rho_m, est_sigma_rho_m))

max <- c(10, #mu alpha 
         12, #delta alpha
         1, #mu rho
         2, #delta rho
         max(true_sigma_alpha_s, est_sigma_alpha_s), max(true_sigma_alpha_m, est_sigma_alpha_m), 
         100, 100) #sigma rho

for (i in 1:length(param_names)){
  params_i <- names[grepl(param_names[i], names, fixed = TRUE) ==T]
  true <- params_i[grepl("true", params_i, fixed = T)]
  est <- params_i[grepl("est", params_i, fixed = T)]
  p <- recov_plot(parameter = param_names[i], true_params = param_df[[true]], 
                  est_params = param_df[[est]], min_val = min[i], max_val = max[i], 
                  plot_name =  recov_titles[i])
  recov_plots[[i]] <- p
}
multi.page <- ggarrange(recov_plots[[1]], recov_plots[[2]], recov_plots[[3]], 
                        recov_plots[[4]], recov_plots[[5]], recov_plots[[6]], 
                        recov_plots[[7]], recov_plots[[8]], nrow=2, ncol=1) 

pdf(file = plot_path) #open pdf file
multi.page ## est and true plots 
param_hist(param_df) #histograms
dev.off()
