# This simulates common normal random variables for use in 
# computing quantiles, power, and sample size. 
library("MASS")

simulateNormal <- function(n_EDTR, n_sim = 1000){
  #Simulates 500 sets of n_sim multivariate normal random variables.
  #
  # Args: 
  #  n_EDTR: number of embedded dynamic treatment regimens (EDTRs).
  #  n_sim: number of normal draws for each of the 500 sets.
  #
  # Returns:
  #  A list of 500 sets, each with n_sim normal random draws.
  mc_list <- vector("list", length = 500)
  for (i in 1:500) {
    mc_list[[i]] <- mvrnorm(n = n_sim, mu = rep(0,n_EDTR), Sigma = diag(n_EDTR))
  }
  return(mc_list)
}
