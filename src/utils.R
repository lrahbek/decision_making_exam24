## ---- Functions for Plotting etc. ---- 

## ---- Recovery Plot: Correlations Between True and Estimated Parameters ----
recov_plot <- function(parameter, true_params, est_params, min_val, max_val, 
                       plot_name){
  
  df <- data.frame(true = true_params, est = est_params)
  RedL <- data.frame(slope = 1, intercept = 0, Correlation='Ideal (slope = 1)')
  regL <- lm(est_params ~ true_params)$coef
  regdf <- data.frame(slope = regL[2], intercept = regL[1],
                      Correlation=sprintf('Actual (slope = %s)', round(regL[2], digits = 3)))
  
  p <- ggplot(df, aes(x = true_params, y = est_params))+
    geom_point(shape = 21)+ xlim(min_val, max_val)+ ylim(min_val, max_val)+
    labs(title = plot_name)+xlab(sprintf("true_%s", parameter))+ ylab(sprintf("est_%s", parameter))+
    geom_abline(aes(slope = slope, intercept = intercept, linetype = Correlation), 
                data = regdf , color = "black")+
    geom_abline(aes(slope = slope, intercept = intercept, linetype = Correlation), 
                data = RedL , color = "red")+
    theme_minimal()
  return(p)
}

## ---- Parameter Histograms: Plots Histograms of Given Parameters ----
param_hist <- function(parameters){
  plot_hist <- list()
  par(mfrow = c(2,2))
  for (i in 1:ncol(parameters)){
    p <- hist(parameters[[i]], main = paste("Histogram of", colnames(parameters)[i]), 
              xlab = colnames(parameters)[i])
    plot_hist[[i]] <- p
  }
  return(plot_hist)
}

## ---- Estimation Error: Correlation Plots Between True-Est Parameters and Vars ----
estimation_error_plot <- function(params, true_params, est_params){
  par(mfrow = c(3, 1))
  sims <- list("contributions" = c_sim, "preferences" = P_sim, "group-belief" = Gb_sim)
  for (i in 1:length(params)){
    parameter <- params[i]
    for (s in 1:length(sims)){
      plot(rowMeans(sims[[s]]), true_params[[i]]-est_params[[i]], 
           xlab = sprintf("Mean %s", names(sims)[s]), 
           ylab = sprintf("true_%1$s - est_%1$s", parameter))
      abline(a = 0, b = 0, col = "red", lwd = 1)
    }
  }
}

## ---- MPD: MAximum Posterior Density Fucntion ----
MPD <- function(x) {
  density(x)$x[which(density(x)$y==max(density(x)$y))]
}