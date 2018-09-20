getOneC <- function(V, Z, alpha = 0.05){
  # Computes a single vector of the ci's given the covariance, a single 
  # vector of standard normal random draws, and type I error rate alpha.
  #
  # Args:
  #  V: covariance matrix of \sqrt{n}(\hat{\theta}).
  #  Z: vector of standard normal random draws.
  #  alpha: type I error rate.
  #
  # Returns:
  #  A single estimate of a nrow(V) length vector the ci's.
  c_alpha <- rep(0,nrow(V))
  
  for (i in 1:nrow(V)) {
    temp <- matrix(0,nrow(Z),nrow(V))
    for (j in setdiff(1:nrow(V),c(i))) {
      temp[,j] <- (Z[,j]-Z[,i])/sqrt(V[i,i]+V[j,j]-2*V[i,j])
    }
    temp <- temp[,-i]
    #temp<-apply(temp,1,sort)
    temp <- apply(temp,1,max)
    c_alpha[i] <- quantile(temp,1-alpha)
  } 
  c_alpha
}

computeC <- function(V, mc_list_SVD, alpha = 0.05) {
  # Computes an estimate for ci (quantiles) averaged 500 times.
  #
  # Args:
  #  V: Covariance matrix of \sqrt{n}(\hat{\theta}).
  #  mc_list_SVD: A list of 500 sets of standard normal random variables.
  #  alpha: Type I Error Probability of excluding optimal EDTR from set of best.
  #
  # Depends:
  #  getOneC(V, Z, alpha).
  #
  # Returns:
  #  An nrow(V) length vector of c'i (quantiles) averaged 500 times.
  c_alpha <- rep(0,nrow(V))
  for (i in 1:500) {
    c_alpha <- c_alpha + getOneC(V, Z = mc_list_SVD[[i]], alpha)
  }
  c_alpha <- c_alpha/500
  c_alpha
}
