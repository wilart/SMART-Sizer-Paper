#Computes the necessary sample size to achieve desired power
# Depends on (1) simulate-normal.R, (2) compute-SVD.R, (3) compute-c-quantiles.R

getOneSampleSize <- function(V, Delta, min_Delta, desired_power, c_alpha, Z){
  # Computes a single estimate for the required sample size to achieve 
  # a specified power.
  #
  # Args:
  #  V: the covariance matrix
  #  Delta: a vector of effect sizes (theta_best-theta_i))_i. The first 0 indicates best DTR
  #  min_Delta: the minimum desired detectable effect size 
  #  desired_power: the desired power
  #  c_alpha: the vector of equicoordinate quantiles for controlling the type I error rate
  #  Z: a simulated normal random variable data set N(0,V)
  #
  # Returns: 
  #  The sample size to achieve desired power

      sample_size <- 0
      N <- min(which(Delta == 0))
      temp <- matrix(0, nrow(Z), nrow(V))
      for (i in which(Delta >= min_Delta)) {
        temp[,i] <- ((Z[,i] - Z[,N]) + c_alpha[i] * sqrt(V[i,i] + V[N,N] - 2 * V[i,N])) / Delta[i]
      }
      temp <- temp[,-N]
      temp <- apply(temp, 1, max)
      sample_size <- (quantile(temp, desired_power))^2
    
  return(sample_size[[1]])
}


computeSampleSize <- function(V, Delta, min_Delta, desired_power, c_alpha, Z = mc_list_SVD) {
  # Averages estimates of necessary sample size over 500 simulated normal data sets
  # see getOnesample_size for more details
  sample_size <- 0
  for (i in 1:500) {
    sample_size <- sample_size + getOneSampleSize(V, Delta, min_Delta, desired_power, c_alpha, Z = Z[[i]])
  }
  sample_size <- sample_size / 500
  return(sample_size)
}
