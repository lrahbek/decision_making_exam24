## ---- PARAMETER RECOVERY: MODEL AR (ALPHA & RHO) ----
install.packages("pacman")
set.seed(42)
pacman::p_load(R2jags, parallel, polspline, tidyverse,ggplot2, extraDistr, ggpubr, truncnorm)
setwd("/work/LauraGivskovRahbek#6659/decision_making_exam24/src/model_AR")
source("../sim_functions.R")
source("../utils.R")

## ---- Simulation Specifications ----
ntokens <- 10 # represents the 10,000
vals <- seq(0, ntokens, 1) # values possible to contribute
groupsize <- 6 # number of participants per group
ngroups <- 4 # for one treatment (e.g. sunni in same-sect groups)
nsubs <- as.integer(groupsize*ngroups)

niter <- 1000

## ---- Empty Arrays for Simulated Data and True and Estimated Parameters ----
true_mu_alpha <- array(NA, c(niter))
true_mu_rho <- array(NA, c(niter))
true_sigma_alpha <- array(NA, c(niter))
true_sigma_rho <- array(NA, c(niter))

est_mu_alpha <- array(NA, c(niter))
est_mu_rho <- array(NA, c(niter))
est_sigma_alpha <- array(NA, c(niter))
est_sigma_rho <- array(NA, c(niter))

c_sim <- array(NA, c(niter, nsubs))
P_sim <- array(NA, c(niter, nsubs))
Gb_sim <- array(NA, c(niter, nsubs))
alpha_sim <- array(NA, c(niter, nsubs))
rho_sim <- array(NA, c(niter, nsubs))

## ---- Run n iterations of Recovery ----
for (i in 1:niter){
  
  # Sample mu and sigma for ALPHA
  mu_alpha <- runif(1, min = 0.001, max = 10)
  sigma_alpha <- rtruncnorm(1, a=0.001, b=Inf, mean = 1, sd = 2)

  # Sample mu and sigma for RHO
  mu_rho <- runif(1, min = 0.001, max = 0.999) #avoid extreme values (0 and 1)
  sigma_rho <- runif(1, min = 10, max = 100)
  
  # Simulate data based on parameters
  sim <- PGGsimAR(nsubs, mu_alpha, mu_rho, sigma_alpha, sigma_rho)
  c <- sim$c
  
  # Save all simulations vars for inspection
  c_sim[i, ] <- c
  P_sim[i, ] <- sim$P
  Gb_sim[i, ] <- sim$Gb
  alpha_sim[i, ] <- sim$alpha
  rho_sim[i, ] <- sim$rho
  
  # Prep and run JAGS model
  data <- list("nsubs", "c") 
  params <- c("mu_alpha", "mu_rho", "sigma_alpha", "sigma_rho")
  
  samples <- jags.parallel(data = data, inits=NULL, parameters.to.save = params, 
                           model.file ="PGG_AR.txt", n.chains=3, n.iter=5000, 
                           n.burnin=1000, n.thin=1, jags.seed=3)
  # Save true and estimated parameters
  true_mu_alpha[i] <- mu_alpha
  true_mu_rho[i] <- mu_rho
  true_sigma_alpha[i] <- sigma_alpha
  true_sigma_rho[i] <- sigma_rho
  
  est_mu_alpha[i] <- MPD(samples$BUGSoutput$sims.list$mu_alpha)
  est_mu_rho[i] <- MPD(samples$BUGSoutput$sims.list$mu_rho)
  est_sigma_alpha[i] <- MPD(samples$BUGSoutput$sims.list$sigma_alpha)
  est_sigma_rho[i] <- MPD(samples$BUGSoutput$sims.list$sigma_rho)
  print(i)
}

## ---- Define vars for saving Data and Plotting Recovery ----
rdata_file <- "../../out/RData/AR_recov_data1000i.RData"
plot_path <- "../../out/AR_recovery_plot1000i.pdf"

param_names <- c("mu_alpha", "mu_rho", "sigma_alpha", "sigma_rho") #names of parameters modeled

recov_titles <- c("Recovery of alpha mode (mu)",  "Recovery of rho mean (mu)", 
                  "Recovery of alpha standard deviation (sigma)", 
                  "Recovery of rho concentration (sigma)")

## ---- Save Data ----
save(true_mu_alpha, true_mu_rho, true_sigma_alpha, true_sigma_rho, 
     est_mu_alpha, est_mu_rho, est_sigma_alpha, est_sigma_rho, 
     c_sim, P_sim, Gb_sim, alpha_sim, rho_sim,
     file = rdata_file)

## ---- Create Plots ----

names <- load(rdata_file, verbose = T) #names of parameter arrays 
param_df <- data.frame("true_mu_alpha" = true_mu_alpha, "est_mu_alpha" = est_mu_alpha, 
                        "true_mu_rho" = true_mu_rho, "est_mu_rho" = est_mu_rho, 
                        "true_sigma_alpha" = true_sigma_alpha, "est_sigma_alpha" = est_sigma_alpha,
                        "true_sigma_rho" = true_sigma_rho, "est_sigma_rho" = est_sigma_rho)

recov_plots <- list() #list to save recovery plots into
#min and max values for plots
min <- c(0, 0, min(true_sigma_alpha, est_sigma_alpha), 
         min(true_sigma_rho, est_sigma_rho))
max <- c(10, 1, max(true_sigma_alpha, est_sigma_alpha), 100)

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
                        recov_plots[[4]], nrow=2, ncol=1) 

pdf(file = plot_path) #open pdf file
multi.page ## est and true plots 
param_hist(param_df) #histograms
estimation_error_plot(params = c("mu_alpha", "mu_rho"), 
                      true_params = list(true_mu_alpha, true_mu_rho), 
                      est_params = list(est_mu_alpha, est_mu_rho))
par(mfrow = c(2,2))
plot(est_mu_rho, true_mu_alpha-est_mu_alpha)
abline(a = 0, b = 0, col = "red", lwd = 1)
plot(true_mu_rho, true_mu_alpha-est_mu_alpha)
abline(a = 0, b = 0, col = "red", lwd = 1)
plot(est_mu_alpha, true_mu_rho-est_mu_rho)
abline(a = 0, b = 0, col = "red", lwd = 1)
plot(true_mu_alpha, true_mu_rho-est_mu_rho)
abline(a = 0, b = 0, col = "red", lwd = 1)

dev.off()
