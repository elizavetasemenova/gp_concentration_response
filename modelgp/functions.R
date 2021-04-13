#-------------------------------------
# kriging
#-------------------------------------

CalcSigma <- function(x1, x2, N_s, eta_, rho_, sigma_, theta_){
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
        if ((x1[j] > theta_[i]) & (x2[k] > theta_[i])){
          Sigma[i, j, k] = Sigma[i, j, k] + eta_[i] ^ 2 * (x1[j] - theta_[i])^2 * (x2[k] - theta_[i])^2 * exp(
            -(x1[j] - x2[k])^ 2 / rho_[i]^ 2)
          
          
        }
      }
    }
  }
  return(Sigma)
}

#-------------------------------------
# inverse logit
#-------------------------------------
inv_logit <- function(x){
  1/ (1 + exp(-x))
}

#-------------------------------------
# mode
#-------------------------------------
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

#-------------------------------------
# DC50
#-------------------------------------

DC50_func <- function(s0, gp, x, eps=0.1){
  
  DC50_1_y <- min(gp) + 0.5* (s0 - min(gp))
  
  diff <- abs(gp - DC50_1_y)
  
  DC50_1 <- min(x[ which(diff - min(diff) < eps)])

  return(list(DC50_1))
}
