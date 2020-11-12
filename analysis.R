  rm(list=ls())
  
  #-------------------------------------
  # Libraries
  #-------------------------------------
  library(tidyverse)
  # library(xtable)
  library(rstan)
  options(mc.cores = parallel::detectCores())
  rstan_options(auto_write = TRUE)
  source("functions.R")
  
  #-------------------------------------
  # Read and prepare data
  #-------------------------------------
  controls <- read.csv(file = 'controls.csv')
  controls <- controls %>% select('Normalized')
  
  compound <- read.csv(file = 'compound.csv')

  compound <- compound %>%
    select(concentration, activity) %>%
    mutate(log_concentration = log10(concentration * 10^(-6))) %>% 
    arrange(log_concentration)
  
  # scaling controls
  sd_global <- sd(compound $activity)
  controls_sc <- controls / sd_global
    
  # plot raw data
  ggplot(compound, aes(x=log_concentration, y=activity)) +
    geom_point() + labs(x ="log concentration (x)", y = "Degradation (y)")

  # add replicate number and unique concentration index
  compound <- compound %>%
    mutate(rep_ind = rep(c(1,2), nrow(compound)/2),
           coord_ind = rep(1: (nrow(compound)/2), each=2),
           activity_sc = activity/ sd_global)
  
  # prepare data for Stan
  x_long <- compound$log_concentration
  x <- unique(x_long)
  N <- length(x)
  coord_ind <- compound$coord_ind
  y_sc <- compound$activity_sc
  N_y <- length(y_sc)
  controls_sc <- controls_sc$Normalized

  stan_dat <- list(N_conc = N, 
                   N_treatment = N_y, 
                   x = x, 
                   y_treatment = y_sc, 
                   x_index = coord_ind,
                   N_control = length(controls_sc),
                   y_control = controls_sc)
  
  #-------------------------------------
  # Fit Stan model
  #-------------------------------------
  
  # compile Stan model
  mod <- stan_model(file="model_invlogit.stan")
  
  niter = 2000
  
  # fit the model
  m <- sampling(mod, 
                data = stan_dat, 
                iter = niter, chains = 4, cores=4)
  
  times <- print(get_elapsed_time(m))
  print(paste0("Stan mean computation time : ", round(mean(apply(times, 1, sum)),2), " sec"))
  print(paste0("Stan sd computation time : ", round(sd(apply(times, 1, sum)),2), " sec"))

  # evaluate convergence by Rhat and n_eff
  print(m, pars = c("sigma_rep",  "sigma", "mu", "rho", "eta"),
        probs = c(0.10, 0.5, 0.75, 0.90),
        digits_summary = 2, include = TRUE)
  
  s <- summary(m, pars = c("sigma_rep",  "sigma", "mu", "rho", "eta"),
               probs = c(0.10, 0.5, 0.90),
               digits_summary = 2, include = TRUE)
  #xtable(s$summary)
  
  # inspect convergence visually
  stan_dens(m, pars = c("sigma_rep",  "sigma", "mu", "rho", "eta"), separate_chains = TRUE)

  #-------------------------------------
  # Kriging
  #-------------------------------------
  
  # extract samples
  mu <- rstan::extract(m, pars="mu")[[1]]
  rho <- rstan::extract(m, pars="rho")[[1]]
  theta <- rstan::extract(m, pars="theta")[[1]]
  sigma <- rstan::extract(m, pars="sigma")[[1]]
  sigma_rep <- rstan::extract(m, pars="sigma_rep")[[1]]
  eta <- rstan::extract(m, pars="eta")[[1]]
  y <- t(rstan::extract(m, pars="y")[[1]])
  N_samp <- length(mu)
  
  N_pred <- 200
  x_pred <- seq(min(x), max(x), length.out=N_pred)
  
  SigmaMM <- CalcSigma(x, x, N_samp, eta, rho, sigma, theta)
  SigmaNM <- CalcSigma(x_pred, x, N_samp, eta, rho, sigma, theta)
  
  SigmaMM_inv <- array(data = NA, dim = c(N_samp, N, N))
  for (k in 1:N_samp){
    SigmaMM_inv[k, , ] = solve(SigmaMM[k, , ])
  }
  
  mean_response <- array(data=NA, dim = c(N_samp, N_pred))
  
  for (k in 1:N_samp){
    mean_response[k,] = mu[k] + ((SigmaNM[k,,] %*% SigmaMM_inv[k,,])  %*% (y[, k] - mu[k]))
  }
  
  theta_mean <- mean(theta)
  
  # convert back to original scale
  gp <- mean_response * sd_global
  gp_mean <- apply(gp, 2, mean)
  gp_95 <- apply(gp, 2, function(x) quantile(x, probs = c(0.025, 0.975)))
  gp.df = data.frame(x = c(x_pred), gp = gp_mean, gp_025 = gp_95[1,], gp_975 = gp_95[2,])
  gp.df = gp.df[order(gp.df$x),]
  
  #-------------------------------------
  # Biologically meaningful parameters
  #-------------------------------------
  
  DC50_samps <- rep(NA, dim(gp)[1])
  Dmax_samps <- rep(NA, dim(gp)[1])
  
  for (i in 1: dim(gp)[1]){
    DC50_samps[i] <- DC50_func(mu[i], gp[i,], gp.df$x, eps=0.15)[[1]]
    Dmax_samps[i] <- min(gp[i,])
  }
  
  DC50 <- quantile(DC50_samps, probs= c(0.025, 0.25, 0.5, 0.75, 0.975))
  Dmax <- quantile(Dmax_samps, probs= c(0.025, 0.25, 0.5, 0.75, 0.975))
  
  h <- hist(DC50_samps, breaks = 500, plot = FALSE)
  px <- h$density
  DC50_draws <- sample(h$mids, size = niter * 1.5, replace = TRUE, prob = px)

  #-------------------------------------
  # Plot
  #-------------------------------------
  
  pdf("figire_one_compound.pdf", height=10, width=18)
  
  par(mai=c(1,1,1,1))
  
  col3 <- "#00998a" 
  plot(0, 
        ylab = 'Target protein degradation (%)', 
        xlab = expression('log'[10]*'-concentration [M]'), 
        pch='',
        ylim= c(min(gp, y), max(gp, y)), 
        xlim= c(min(x, x_pred), max(x, x_pred)),
        cex.lab = 2.0)
   
  lims <- par("usr")

  # controls 
  points(rep(min(x_long), length(controls_sc)) + rnorm(length(controls_sc), 0, 0.1), controls$Normalized, pch=19, col='grey40')
  
  # observations
  points(x_long, y_sc * sd_global, pch=15, col="tomato", cex=1.5)
  
  # theta
  abline(v=theta_mean, col="coral3", lty = 2, lwd=2) # coral3
  theta_95 <- quantile(theta, probs = c(0.025, 0.25, 0.75, 0.975))
  rect(theta_95[1], lims[3], theta_95[4], lims[4],  col=alpha("coral3", 0.3), border=FALSE)
  text(theta_mean+0.5, lims[3]+5, expression('log'[10]*'(PoD)'), cex=2.0)
  
  # Dmax
  rect(lims[1], Dmax[1], lims[2], Dmax[5],  col=alpha("goldenrod1", 0.4), border=FALSE)
  abline(h=Dmax[3], col="goldenrod4", lty = 2, lwd=2) # darkgreen
  text(lims[1]+0.3, Dmax[3]+5,  expression('D'["max"]), cex=2.0)
  
  # DC50
  rect(DC50[1], lims[3], DC50[5], lims[4],  col=alpha("steelblue", 0.3), border=FALSE)
  abline(v=DC50[3], col="steelblue", lty= 2, lwd=2)
  text(DC50[3]+0.5, lims[3]+5, expression('log'[10]*'(DC'[50]*')'), cex=2.0)
  
  polygon(c(gp.df$x,rev(gp.df$x)),c(gp.df$gp_975,rev(gp.df$gp_025)),
          col=alpha(col3,0.3), border=NA)

  lines(gp.df$x, gp.df$gp, pch=19, col="darkslategray", lwd=2)

  legend("topright", 
          legend = c("observations", 
                     "controls",
                     #"GP draws", 
                     "GP mean",
                     "GP 95% BCI", 
                     expression('log'[10]*'(PoD)'), 
                     expression('log'[10]*'(PoD) 95% BCI'),
                     expression('D'["max"]),
                     expression('D'["max"]*' 95% BCI'),
                     expression('log'[10]*'(DC'[50]*')'),
                     expression('log'[10]*'(DC'[50]*')'*' 95% BCI')),
          col = c("tomato",
                  "grey40",
                  "darkslategray", 
                  alpha(col3,0.8), 
                  "coral3", 
                  "coral3", 
                  "goldenrod1",
                  "goldenrod1",
                  "steelblue",
                  alpha("steelblue", 0.6)),
         pch = c(15,19, NA, 15, NA, 15, NA, 15, NA, 15), 
         lty=c(NA,NA, 1,NA,2,NA, 2, NA, 2, NA),
         lwd=c(NA,NA,3,NA,3,3,3, NA, 3, NA),
         bty = "n", 
         pt.cex = 1, 
         cex = 1.4, 
         text.col = "black",
         inset = c(0.05, 0.01))
  
  dev.off()

  