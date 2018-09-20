#Estimate the true EDTR outcomes
Design1Simulate <- function(n, delta, n_sim){
  # Source: ADHD simulation code from Prof. Ashkan Ertefaie
  # Simulates n_sim design 1 SMART trials with n individuals with delta difference between best and non-best (confirm)
  #
  # Args: 
  #  n: Number of individuals in each SMART study simulation
  #  delta: Effect size
  #  n_sim: Number of simulations
  #
  # Returns: 
  # A list of n_sim SMART trial simulations each consisting of a data frame of the simulated variables, and outcomes
  ADHD.sim.list <- vector("list", length = n_sim)
  for (i in 1:n_sim){
    o11 <- rnorm(n,0,1)
    o12 <- rnorm(n,0,1)
    a1 <- 2*rbinom(n,1,.5)-1 
    o21<-rnorm(n,0.5*o11+0.5*I(a1==-1),1)
    o22<-rnorm(n,0.5*o12+0.5*I(a1==+1),1)
    r<- (o21>0) #Indicator of Response
    
    s<-rep(0,n)
    s[r==0]<-1
    
    a2<-2*rbinom(n,1,0.5)-1
    
    #####################
    y1 <- (1+o11-o12+o22+o21+a1*(delta+.5*o11) +s*(a2*delta/2 ))
    
    y <- y1 + rnorm(n,0,1)
    
    ADHD.sim.list[[i]]<-data.frame(o11,o12,a1,o21,o22,r,s,a2,y)
  }
  ADHD.sim.list
}

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
  return(theta.hat)
}

set.seed(1237)
sim.asymptotic <- Design1Simulate(n = 10000, delta = 0.25, n_sim = 1000)
library(doParallel)
no_cores <- detectCores()
cl <- makeCluster(no_cores)
registerDoParallel(cl)
clusterExport(cl,"Design1ComputeBeta",envir=environment())

beta.true.list <- parLapply(cl, sim.asymptotic, function(x) t(Design1ComputeBeta(dat = x, n = 10000)))
stopCluster(cl)
beta.true <- colMeans(do.call(rbind,beta.true.list))
round(beta.true,3)
theta.true <- round(Design1ComputeTheta(beta.true),3)
theta.true
Delta.DR <- max(theta.true)-theta.true

