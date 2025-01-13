## ---- PARAMETER ESTIMATION: MODEL A (ALPHA) ----
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
SeTr <- unique(pgg_data$SeTr) #sect and treatment
class <- unique(pgg_data$class) #class
nsubs_ls <- list()
for (i in 1:length(SeTr)){
  nsubs_ls[i] <- nrow(pgg_data[pgg_data$SeTr == SeTr[i],])
  names(nsubs_ls)[i] <- SeTr[i]
}

file_names <- list()

## ---- Run Hierarchical Parameter Estimation on each Sect & Treatment ----
for (i in 1:length(SeTr)){
  
  c <- (pgg_data$c[pgg_data$SeTr == SeTr[i]])/1000 #rescale to be between 0 - 10
  nsubs <- nsubs_ls[[i]]
  data <- list("nsubs", "c") 
  params <- c("mu_alpha", "sigma_alpha")
  
  samples <- jags.parallel(data = data, inits=NULL, parameters.to.save = params, 
                           model.file ="PGG_A.txt", n.chains=3, n.iter=5000, 
                           n.burnin=1000, n.thin=1, jags.seed=3)
  
  file_name <- sprintf("../../out/EstData/param_est_A_%s.RData", SeTr[i])
  save(samples, file = file_name)
  print(sprintf("Samples saved at %s", file_name))
  file_names <- append(file_names, file_name)
}

save(SeTr, file_names, file = "../../out/EstData/est_vars_A.RData")


