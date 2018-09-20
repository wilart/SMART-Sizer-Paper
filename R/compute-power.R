#Defines functions for computing the power using Monte Carlo simulation
#Depends on (1) simulate-normal and (2) compute-SVD to generate N(0,V) random variables
#Depends on (3) compute-c-quantiles to compute the equicoordinate quantiles 
#which control the type I error rate
getOnePower <- function(V, Delta, min_Delta, sample_size, c_alpha, Z){
  # Computes the power using a data set of simulated normal random variables N(0,V)
  #
  # Args:
  #  V: the covariance matrix
  #  Delta: a vector of effect sizes (theta_best-theta_i))_i
  #  min_Delta: the minimum desired detectable effect size with first 0 indicating best EDTR.
  #  sample_size: the sample size
  #  c_alpha: the vector of equicoordinate quantiles for controlling the type I error rate
  #  Z: a simulated normal random variable data set N(0,V)

  power <- 0
  N <- min(which(Delta == 0))
  temp <- matrix(0, nrow(Z), nrow(V))
  for (i in which(Delta >= min_Delta)) {
    temp[,i] <- ((Z[,i] - Z[,N]) + c_alpha[i] * sqrt(V[i,i] + V[N,N] - 2 * V[i,N])) / Delta[i]
  }
  temp <- temp[,-N]
  temp <- apply(temp, 1, max)
  power <- mean(temp < sqrt(sample_size))
    
  return(power)
}

computePower <- function(V, Delta, min_Delta, sample_size, c_alpha, Z = mc_list_SVD) {
  ##Averages the power over 500 normal random variable data sets
  ##See documentation for getOnePower for more information
  ## Z is a list of N(0,V) simulated random variables matrices
  avg_power <- 0
  for (i in 1:500) {
    avg_power <- avg_power + getOnePower(V, Delta, min_Delta, sample_size, c_alpha, Z = Z[[i]])
  }
  avg_power <- avg_power / 500
  return(avg_power)
}
