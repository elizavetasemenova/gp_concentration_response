  rm(list=ls())
  
  # set working directory 
  #setwd(path_to_directory)
  
  #-------------------------------------
  # Include libraries
  #-------------------------------------
  library(utils)
  library(tidyverse)
  library(ggplot2)
  library(rstan)
  options(mc.cores = parallel::detectCores())
  rstan_options(auto_write = TRUE)
  
  
  source("modelgp/functions.R")
  
  # number of iteration of an MCMC (niter) and the number of prediction points (N_pred)
  # define computation time
  niter <- 10000
  N_pred <- 200
  
  #-------------------------------------
  # Read and plot data
  #-------------------------------------
  controls <- read.csv(file = 'datagp/controls.csv')
  controls <- controls %>% select('Normalized')
  
  compound <- read.csv(file = 'datagp/compound.csv')
  
  head(controls)
  
  head(compound)
  

  compound <- compound %>%
    select(concentration, activity) %>%
    mutate(log_concentration = log10(concentration * 10^(-6))) %>% 
    arrange(log_concentration)
  
  # scaling controls
  sd_global <- sd(compound $activity)
  controls_sc <- controls / sd_global
    
  # plot raw data
  # create an artificial dataset to plot controls, and add jitter to their x-coordinate
  controls_to_plot <- controls
  controls_to_plot['x'] = rep(min(compound$log_concentration) , nrow(controls_to_plot)) + rnorm(nrow(controls_to_plot), 0,0.08)
  names(controls_to_plot) <- c("activity", "log_concentration")
  
  p <-  ggplot(compound, aes(x=log_concentration, y=activity)) +
    geom_point(colour="tomato") + labs(x ="log concentration (x)", y = "Degradation (y)") +
    geom_point(data=controls_to_plot, color="grey40")
  
  p
 
  
  #-------------------------------------
  # Prepare data for Stan
  #-------------------------------------
  
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
  # Compile Stan model
  #-------------------------------------
  
  mod <- stan_model(file="modelgp/model_invlogit_nu_g.stan")
  
  #-------------------------------------
  # Fit Stan model
  #-------------------------------------
  
  # fit the model
  m <- sampling(mod, 
                data = stan_dat, 
                iter = niter, chains = 4, cores=4)
  
  # to avoid running all of the below code, call
  # post_GP(m)
  
  times <- print(get_elapsed_time(m))
  print(paste0("Stan mean computation time : ", round(mean(apply(times, 1, sum)),2), " sec"))
  print(paste0("Stan sd computation time : ", round(sd(apply(times, 1, sum)),2), " sec"))
  
  #-------------------------------------
  # evaluate convergence
  #-------------------------------------
  
  pars = c("sigma_rep",  "sigma", "mu", "rho", "eta", "nu", "g")
  
  # evaluate convergence by Rhat and n_eff
  print(m, 
        pars = pars,
        probs = c(0.10, 0.5, 0.75, 0.90),
        digits_summary = 2, include = TRUE)
  
  # s <- summary(m, pars = c("sigma_rep",  "sigma", "mu", "rho", "eta"),
  #              probs = c(0.10, 0.5, 0.90),
  #              digits_summary = 2, include = TRUE)
  #xtable(s$summary)
  
  # traceplots
  traceplot(m, pars = pars, inc_warmup=FALSE)
  
  # posterior density plots
  stan_dens(m, pars = pars, separate_chains = TRUE)
  
  #-------------------------------------
  # Kriging
  #-------------------------------------
  
  mu <- rstan::extract(m, pars="mu")[[1]]
  rho <- rstan::extract(m, pars="rho")[[1]]
  theta <- rstan::extract(m, pars="theta")[[1]]
  sigma <- rstan::extract(m, pars="sigma")[[1]]
  sigma_rep <- rstan::extract(m, pars="sigma_rep")[[1]]
  eta <- rstan::extract(m, pars="eta")[[1]]
  g_par <- rstan::extract(m, pars="g")[[1]]
  y <- t(rstan::extract(m, pars="y")[[1]])
  N_samp <- length(mu)
  
  x_pred <- seq(min(x), max(x), length.out=N_pred)
  
  SigmaMM <- CalcSigma_g(x, x, N_samp, eta, rho, sigma, theta, g_par)
  SigmaNM <- CalcSigma_g(x_pred, x, N_samp, eta, rho, sigma, theta, g_par)
  
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
  
  #-------------------------------------
  # Plot results
  #-------------------------------------

  # Uncomment this line to print the final plot into an external pdf.
  #pdf("figure_one_compound.pdf", height=10, width=18)
  
  par(mai=c(1,1,1,1))
  ptscex <- 2
  title <- ""
  
  col3 <- "#00998a" 
  #point_cols <- c("tomato", "dodgerblue4", "darkslategray", "goldenrod3")
  point_cols <- c("tomato", "dodgerblue4", "darkolivegreen3", "goldenrod3")
  
  plot(0, 
       ylab = 'Target protein degradation (%)', 
       xlab = expression('log'[10]*'-concentration [M]'), 
       pch='',
       #ylim= c(min(gp, y), max(gp, y)), 
       ylim=c(-110, 10),
       xlim= c(min(x, x_pred), max(x, x_pred)),
       cex.lab = 2.0,
       main=title,
       cex.main=2)
  
  lims <- par("usr")
  
  # controls 
  points(rep(min(x_long), length(controls_sc)) + rnorm(length(controls_sc), 0, 0.1), controls$Normalized, pch=8, col='grey40', cex=ptscex)
  
  # observations
  #points(x_long, y_sc * sd_global, pch=15, col="tomato", cex=1.5)
  points(x_long, y_sc * sd_global, pch=19, col=point_cols[compound$rep_ind], cex=ptscex)
  
  # theta
  abline(v=theta_mean, col="coral3", lty = 2, lwd=2) # coral3
  theta_95 <- quantile(theta, probs = c(0.025, 0.25, 0.75, 0.975))
  rect(theta_95[1], lims[3], theta_95[4], lims[4],  col=alpha("coral3", 0.2), border=FALSE)
  text(theta_mean+0.5, lims[3]+5, expression('log'[10]*'(PoD)'), cex=2.0)
  
  
  # Dmax
  rect(lims[1], Dmax[1], lims[2], Dmax[5],  col=alpha("goldenrod1", 0.4), border=FALSE)
  abline(h=Dmax[3], col="goldenrod4", lty = 2, lwd=2) # darkgreen
  text(lims[1]+0.3, Dmax[3]+5,  expression('D'["max"]), cex=2.0)
  
  
  # DC50
  rect(DC50[1], lims[3], DC50[5], lims[4],  col=alpha("steelblue", 0.2), border=FALSE)
  abline(v=DC50[3], col="steelblue", lty= 2, lwd=2)
  text(DC50[3]+0.5, lims[3]+5, expression('log'[10]*'(DC'[50]*')'), cex=2.0)
  
  polygon(c(gp.df$x,rev(gp.df$x)),c(gp.df$gp_975,rev(gp.df$gp_025)),
          col=alpha(col3,0.2), border=NA)
  
  lines(gp.df$x, gp.df$gp, pch=19, col="darkslategray", lwd=2)
  
  legend("topright", 
         legend = c("observations-rep1", 
                    "observations-rep2", 
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
         col = c(point_cols[1],
                 point_cols[2],
                 "grey40",
                 "darkslategray", 
                 alpha(col3,0.8), 
                 "coral3", 
                 "coral3", 
                 "goldenrod1",
                 "goldenrod1",
                 "steelblue",
                 alpha("steelblue", 0.6)),
         pch = c(19, 19, 8, NA, 15, NA, 15, NA, 15, NA, 15), 
         lty=c(NA,NA, NA, 1,NA,2,NA, 2, NA, 2, NA),
         lwd=c(NA,NA, NA,3,NA,3,3,3, NA, 3, NA),
         bty = "n", 
         pt.cex = 1, 
         cex = 1.4, 
         text.col = "black",
         inset = c(0.05, 0.01))
  
  # Uncomment this line to print the final plot into an external pdf.
  #dev.off()
  
  #-------------------------------------
  # Print estimates of biological parameters
  #-------------------------------------
  
  out.df <- data.frame(Dmax = c(Dmax[1], Dmax[3], Dmax[5]),
                       DC50 = c(DC50[1], DC50[3], DC50[5]),
                       theta = c(theta_95[1], theta_mean, theta_95[4]))
  
  cat('Parameter estimates :\n')
  round(t(out.df),2)
  
  #-------------------------------------
  # Computational environment
  #-------------------------------------
  sessionInfo()
  
  
