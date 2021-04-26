
#---------------------------------------------------------------------------------------------------------------
# kriging: CalcSigma_g
#---------------------------------------------------------------------------------------------------------------
CalcSigma_g <- function(x1, x2, N_s, eta_, rho_, sigma_, theta_, k_){
  N1 <- length(x1)
  N2 <- length(x2)
  Sigma <- array(data=0, dim = c(N_s, N1, N2))
  sigma2 <- sigma_ * sigma_
  
  if (N1 == N2){
    for (i in 1:N1)
      Sigma[, i, i] = sigma2[i]
  }
  
  for (i in 1:N_s){
    for (j in 1:N1){
      for (k in 1:N2){
        
        d1 = x1[j] - theta_[i]
        d2 = x2[k] - theta_[i]
        K2 = eta_[i]^2 * d1^2 * d2^2 * exp(-(d1-d2)^2/(rho_[i]^2))
        
        Sigma[i, j, k] = Sigma[i, j, k] +  inv_logit(k_[i]*d1)* K2 *inv_logit(k_[i]*d2)
        
      }
    }
  }
  return(Sigma)
}

#---------------------------------------------------------------------------------------------------------------
# post-process GP
#---------------------------------------------------------------------------------------------------------------
post_GP <- function(m, N_pred=200, legend=TRUE, title="", ptscex=1){
  times <- print(get_elapsed_time(m))
  print(paste0("Stan mean computation time : ", round(mean(apply(times, 1, sum)),2), " sec"))
  print(paste0("Stan sd computation time : ", round(sd(apply(times, 1, sum)),2), " sec"))
  
  
  # evaluate convergence
  pars = c("sigma_rep",  "sigma", "mu", "rho", "eta", "nu", "g")
  
  # evaluate convergence by Rhat and n_eff
  print(m, 
        pars = pars,
        probs = c(0.10, 0.5, 0.75, 0.90),
        digits_summary = 2, include = TRUE)
  
  # Kriging
  
  mu <- rstan::extract(m, pars="mu")[[1]]
  rho <- rstan::extract(m, pars="rho")[[1]]
  theta <- rstan::extract(m, pars="theta")[[1]]
  sigma <- rstan::extract(m, pars="sigma")[[1]]
  sigma_rep <- rstan::extract(m, pars="sigma_rep")[[1]]
  eta <- rstan::extract(m, pars="eta")[[1]]
  g_par <- rstan::extract(m, pars="g")[[1]]
  #nu <- rstan::extract(m, pars="nu")[[1]]
  y <- t(rstan::extract(m, pars="y")[[1]])
  N_samp <- length(mu)
  
  #N_pred <- 200
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
  
  
  #par(mai=c(1,1,1,1))
  
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
  
  if(legend==TRUE){
  if (max(compound$rep_ind)==4){
    legend("topright", 
           legend = c("observations-rep1", 
                      "observations-rep2", 
                      "observations-rep3", 
                      "observations-rep4", 
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
                   point_cols[3],
                   point_cols[4],
                   "grey40",
                   "darkslategray", 
                   alpha(col3,0.8), 
                   "coral3", 
                   "coral3", 
                   "goldenrod1",
                   "goldenrod1",
                   "steelblue",
                   alpha("steelblue", 0.6)),
           pch = c(19, 19, 19, 19,8, NA, 15, NA, 15, NA, 15, NA, 15), 
           lty=c(NA,NA, NA, NA,NA, 1,NA,2,NA, 2, NA, 2, NA),
           lwd=c(NA,NA, NA, NA,NA,3,NA,3,3,3, NA, 3, NA),
           bty = "n", 
           pt.cex = 1, 
           cex = 1.4, 
           text.col = "black",
           inset = c(0.05, 0.01))
    
  } 
  
  
  if (max(compound$rep_ind)==3){
    legend("topright", 
           legend = c("observations-rep1", 
                      "observations-rep2", 
                      "observations-rep3", 
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
                   point_cols[3],
                   "grey40",
                   "darkslategray", 
                   alpha(col3,0.8), 
                   "coral3", 
                   "coral3", 
                   "goldenrod1",
                   "goldenrod1",
                   "steelblue",
                   alpha("steelblue", 0.6)),
           pch = c(19, 19, 19,8, NA, 15, NA, 15, NA, 15, NA, 15), 
           lty=c(NA,NA, NA, NA, 1,NA,2,NA, 2, NA, 2, NA),
           lwd=c(NA,NA, NA, NA,3,NA,3,3,3, NA, 3, NA),
           bty = "n", 
           pt.cex = 1, 
           cex = 1.4, 
           text.col = "black",
           inset = c(0.05, 0.01))
    
  } 
  
  if (max(compound$rep_ind)==2){
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
    
  } 
  
  if (max(compound$rep_ind)==1){
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
           pch = c(19,8, NA, 15, NA, 15, NA, 15, NA, 15), 
           lty=c(NA,NA, 1,NA,2,NA, 2, NA, 2, NA),
           lwd=c(NA,NA,3,NA,3,3,3, NA, 3, NA),
           bty = "n", 
           pt.cex = 1, 
           cex = 1.4, 
           text.col = "black",
           inset = c(0.05, 0.01))
    } 
  }
  
  
  # Uncomment this line to print the final plot into an external pdf.
  #dev.off()
  
  #-------------------------------------
  # Print parameter estimates
  #-------------------------------------
  
  out.df <- data.frame(Dmax = c(Dmax[1], Dmax[3], Dmax[5]),
                       DC50 = c(DC50[1], DC50[3], DC50[5]),
                       theta = c(theta_95[1], theta_mean, theta_95[4]))
  
  cat('Parameter estimates :\n')
  round(t(out.df),2)
  
  # library(xtable)
  # xtable(t(out.df))
  
}

#---------------------------------------------------------------------------------------------------------------
# inverse logit
#---------------------------------------------------------------------------------------------------------------
inv_logit <- function(x){
  1/ (1 + exp(-x))
}

#---------------------------------------------------------------------------------------------------------------
# mode
#---------------------------------------------------------------------------------------------------------------
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

#---------------------------------------------------------------------------------------------------------------
# DC50
#---------------------------------------------------------------------------------------------------------------

DC50_func <- function(s0, gp, x, eps=0.1){
  
  DC50_1_y <- min(gp) + 0.5* (s0 - min(gp))
  
  diff <- abs(gp - DC50_1_y)
  
  DC50_1 <- min(x[ which(diff - min(diff) < eps)])

  return(list(DC50_1))
}
