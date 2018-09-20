# This compute the singular value decomposition (SVD) over a list of 
# simulated data sets. 
# We adapted code from the mvrnorm function in the MASS package to deal with empirical 
# covariance matrices outside the parameter space
# (non-positive definite matrices)
# Venables, W. N. & Ripley, B. D. (2002) Modern Applied Statistics with S. Fourth Edition. 
# Springer, New York. ISBN 0-387-95457-0


ComputeSVD <- function(V, mc_list, tol = 1e-6){
  # Returns the transformed normal random variables N(0,I) to N(0,V) using SVD.
  # 
  # Args: 
  #  V: covariance matrix of \sqrt{n}(\hat{\theta}).
  #  mc_list: list of simulated N(0,I) random variables
  #  tol: tolerance for eigenvalues to be negative
  #
  # Returns:
  #  A list of N(0,V) random variables
  eS <- eigen(V, symmetric = TRUE)
  ev <- eS$values
  if (!all(ev >= -tol * abs(ev[1L])))
    stop("'Sigma' is not positive definite") #MVRNORM code
  
  lapply(mc_list, FUN = function(x) t(eS$vectors %*%
                                        (diag(sqrt(pmax(ev,0)))   
                                         %*% t(x)))) #see MVRNORM
}
