# Computes set of best for each data set and empirical power
ComputeBestSet <- function(theta, c_alpha, V, sample_size){
  ##Returns a vector of indicators of membership to the set of best DTRs
  in.best.set.indicator <- rep(1,nrow(V))
  for(i in 1:nrow(V)){
    for(j in 1:nrow(V)){
      in.best.set.indicator[i] <- in.best.set.indicator[i]*(sqrt(sample_size)*theta[i] >= sqrt(sample_size)*theta[j]-c_alpha[i]*sqrt(V[i,i]+V[j,j]-2*V[i,j]))
    }
    
  }  
  return(in.best.set.indicator)
}


##Compute the empirical power
ComputeEmpiricalPower <- function(Delta, theta.true, in.best.set.indicator.one){
  not.best.matrix <- 1 - do.call(rbind,in.best.set.indicator.one)
  true.not.best <- which(round(max(theta.true)-theta.true,2)>=Delta)
  exclude.indicator <- 1
  for (i in 1:length(true.not.best)){
    exclude.indicator <- exclude.indicator*not.best.matrix[,true.not.best[i]]
  }
  return(mean(exclude.indicator))
}

