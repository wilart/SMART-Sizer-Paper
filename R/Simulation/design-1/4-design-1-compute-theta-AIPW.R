#Simulation Design 1: Computes betas and thetas using AIPW for the first SMART design example.
Design1ComputeBeta <- function(dat, n){
  n <- nrow(dat)
  dat$ID<-1:n
  
  dat.R1 <- dat.R2 <- dat[dat$r==1,];
  dat.R1$a2 <- 1; dat.R2$a2 <- -1;
  
  rep.dat <- rbind(dat.R1, dat.R2, dat[dat$r==0,])
  w <- rep(1,dim(rep.dat)[1])
  w[rep.dat$s==1]<-2
  
  y<-rep.dat$y
  
  des <- cbind(1,rep.dat$a1,rep.dat$a2)
  a1 <- dat$a1
  a2 <- dat$a2
  y2 <- dat$y
  
  
  den <- num <- 0
  AI.m<-matrix(c(1,1,-1,1,1,-1,-1,-1),ncol=2,byrow=T)
  
  for(k in 1:4){
    
    a1i <- AI.m[k,1]
    a2j <- AI.m[k,2]
    ####### BELOW IS FOR POINT ESTIMATE
    cof <- coef(lm(y~o11+o12+o21+o22+a1+a1*o11+a2,data= rep.dat,weight=w) )
    
    des <- cbind(1,rep(a1i, n),rep(a2j,n))
    
    phi2 <- cof[1]+cof[2]*dat$o11+cof[3]*dat$o12+ cof[4]*dat$o21+cof[5]*dat$o22+cof[6]*a1i+1*cof[7]*a2j+cof[8]*a1i*dat$o11
    
    phi22 <- cof[1]+cof[2]*rep.dat$o11+cof[3]*rep.dat$o12+ cof[4]*rep.dat$o21+cof[5]*rep.dat$o22+cof[6]*rep.dat$a1+1*cof[7]*a2j+cof[8]*rep.dat$a1*rep.dat$o11
    
    cof2 <- coef(lm(phi22~o11+o12+a1+a1*o11,data= rep.dat,weight=w) )
    phi1 <- cof2[1]+cof2[2]*dat$o11+cof2[3]*dat$o12+cof2[4]*a1i+cof2[5]*a1i*dat$o11
    deni <- t(des)%*%(des)
    numi <- (1)*t(dat$r*des)%*%((y2-phi2+1*cof[7]*a2j)*((a1==a1i)/mean(a1==a1i)))+t(1*(1-dat$r)*des)%*%((y2-phi2)*((a1==a1i&a2==a2j)/mean(a1==a1i&a2==a2j)))+t(dat$r*des)%*%((phi2-phi1-1*cof[7]*a2j)*((a1==a1i)/mean(a1==a1i)))+t((1-dat$r)*des)%*%((phi2-phi1)*((a1==a1i)/mean(a1==a1i)))+t(1*des)%*%(phi1)
    
    den <- den+deni
    num <- num+numi
    #########################END of POINT ESTIMATION
  }
  beta<-solve(den)%*%num # THIS IS THE POINT ESTIMATE USING THE MORE EFFICIENT ESTIMATOR. (\BETA AUGMENTED)
  return(beta)
}
Design1ComputeTheta <- function(beta){
  #(1,1)
  e1 <- t(beta)%*%(c(1,1,1))
  
  #(-1,1)
  e2 <- t(beta)%*%(c(1,-1,1))
  
  #(1,-1)
  e3 <- t(beta)%*%(c(1,1,-1))
  
  #(-1,-1)
  e4 <- t(beta)%*%(c(1,-1,-1))
  
  # best txt regime: (1,1,-1)
  theta.hat <- cbind(e1,e2,e3,e4)
}



n_grid <- seq(100, 500, 50)

design.1.beta.by.n.delta.0.25 <- vector("list", length = length(n_grid))

library(doParallel)
no_cores <- detectCores()
cl <- makeCluster(no_cores)
registerDoParallel(cl)

#Compute the empirical AIPW betas from the simulated data sets
design.1.beta.by.n.delta.0.25 <- foreach(i = 1:length(n_grid)) %:% 
  foreach(j = 1:1000) %dopar% {
    Design1ComputeBeta(design.1.sim.by.n.delta.0.25[[i]][[j]], n_grid[i])
  }

stopCluster(cl)


design.1.theta.by.n.delta.0.25 <- vector("list", length = length(n_grid))
#Transform the betas to thetas
for (i in 1:length(n_grid)){
  design.1.theta.by.n.delta.0.25[[i]] <- lapply(design.1.beta.by.n.delta.0.25[[i]], Design1ComputeTheta)
}
