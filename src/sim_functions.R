## ---- Functions for Simulations ---- 

## ---- Model-AR: Simulates a Group based on MU and SIGMA for ALPHA and RHO ----
PGGsimAR <- function(nsubs, mu_alpha, mu_rho, sigma_alpha, sigma_rho){
  # Define arrays for data
  c <- array(NA, c(nsubs)) #contribution 
  P <- array(NA, c(nsubs)) #preference
  Gb <- array(NA, c(nsubs)) #group belief
  alpha <- array(NA, c(nsubs)) #alpha param for each sub
  rho <- array(NA, c(nsubs)) #rho param for each sub
  
  # Define rate and shape from treatment level mu and sigma for alpha
  rate_a <- (mu_alpha + sqrt(mu_alpha^2 + 4*(sigma_alpha)^2)) / (2*(sigma_alpha)^2)
  shape_a <- 1 + (mu_alpha * rate_a)
  
  # Define shapes from treatment level mu and sigma for rho
  shape_r1 <- mu_rho * sigma_rho
  shape_r2 <- (1 - mu_rho) * sigma_rho
  
  for (s in 1:nsubs){
    # Define subject level parameters
    alpha[s] <- rgamma(1, shape = shape_a, rate = rate_a)
    rho[s] <- rbeta(1, shape_r1, shape_r2)
    
    # Define behaviour in PGG
    Gb[s] <- rpois(1, alpha[s]) #group belief
    P[s] <- rho[s] * Gb[s] #preferences
    c[s] <-  extraDistr::rtpois(1, P[s], a = 0, b = 10) #contribution first trial
    
  }
  
  results <- list(c = c, P = P, Gb = Gb, alpha = alpha, rho = rho)
  return(results)
}

## ---- MODEL-A: Simulates a Group based on MU and SIGMA for ALPHA ----
PGGsimA <- function(nsubs, mu_alpha, sigma_alpha, rho){
  # Define arrays for data
  c <- array(NA, c(nsubs)) #contribution 
  P <- array(NA, c(nsubs)) #preference
  Gb <- array(NA, c(nsubs)) #group belief
  alpha <- array(NA, c(nsubs)) #alpha param for each sub

  # Define rate and shape from treatment level mu and sigma for alpha
  rate_a <- (mu_alpha + sqrt(mu_alpha^2 + 4*(sigma_alpha)^2)) / (2*(sigma_alpha)^2)
  shape_a <- 1 + (mu_alpha * rate_a)
  
  
  for (s in 1:nsubs){
    # Define subject level parameters
    alpha[s] <- rgamma(1, shape = shape_a, rate = rate_a)

    # Define behaviour in PGG
    Gb[s] <- rpois(1, alpha[s])
    P[s] <- rho * Gb[s] #preferences first trial
    c[s] <-  extraDistr::rtpois(1, P[s], a = 0, b = 10) #contribution first trial
  }
  results <- list(c = c, P = P, Gb = Gb, alpha = alpha)
  return(results)
}
