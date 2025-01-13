## ---- PARAMETER COMPARE: MODEL A (ALPHA) ----
install.packages("pacman")
set.seed(42)
pacman::p_load(R2jags, parallel, polspline, tidyverse,ggplot2, extraDistr, ggpubr)
setwd("/work/LauraGivskovRahbek#6659/decision_making_exam24/src/model_A")
source("../utils.R")

## ---- Load cleaned data ----
pgg_data <- read_csv("../../in/data_clean.csv")

## ---- Characteristics ----
sects <- unique(pgg_data$sect) # C, Sh, S
treatments <- unique(pgg_data$treatment) # 0 (same-sect), 1 (mixed-sect)
class <- unique(pgg_data$class) #class

file_names <- list()

## ---- Run Hierarchical Parameter Comparison on each Sect & Treatment ----
for (i in 1:length(sects)){
  
  c_s <- (pgg_data$c[pgg_data$sect == sects[i] & pgg_data$treatment == 0])/1000 #rescale to be between 0 - 10
  c_m <- (pgg_data$c[pgg_data$sect == sects[i] & pgg_data$treatment == 1])/1000 #rescale to be between 0 - 10
  nsubs_s <- length(c_s)
  nsubs_m <- length(c_m)
  
  data <- list("nsubs_s", "nsubs_m", "c_s", "c_m") 
  params <- c("mu_alpha", "delta_alpha", "sigma_alpha_s", "sigma_alpha_m")
  
  samples <- jags.parallel(data = data, inits=NULL, parameters.to.save = params, 
                           model.file ="PGGcompare_A.txt", n.chains=3, n.iter=5000, 
                           n.burnin=1000, n.thin=1, jags.seed=3)
  
  file_name <- sprintf("../../out/EstData/param_est_compare_A_%s.RData", sects[i])
  save(samples, file = file_name)
  print(sprintf("Samples saved at %s", file_name))
  file_names <- append(file_names, file_name)
}

save(sects, file_names, file = "../../out/EstData/est_compare_vars_A.RData")

