#Simulation Design 2: Simulate datasets
delta = 2
n_sim = 1000

##From Autism-IPW-Table2.R
Design2Simulate <- function(n, delta, n_sim){
  design2.sim.list <- vector("list", length = n_sim)
  for (i in 1:n_sim){
    
    o11<-rnorm(n,0,1)
    
    o12<-rnorm(n,0,1)
    
    a1<-2*rbinom(n,1,.5)-1
    
    o21<-rnorm(n,0.5*o11 + 0.5*I(a1==-1),1)
    
    o22<-rnorm(n,0.5*o12 + 0.5*I(a1==+1),1)
    
    r<- (o22>0) 
    
    s<-rep(0,n)
    s[r==0]<-1 # Non-Responder indicator
    
    
    a2<-sample(c(1,2,3,4),n,replace=TRUE)
    
    
    #####################
    y1<- (1+o11-o12+o22+o21+I(a1==-1)*(delta+.5*o11) +1*s*I(a1==-1)*(-delta/4*I(a2==1)+delta/2*I(a2==2)+0*I(a2==3)+delta/2*o21*I(a2==2) )) 

    y<-y1+rnorm(n,0,1)
    design2.sim.list[[i]]<-data.frame(o11,o12,a1,o21,o22,r,s,a2,y)
  }
  
  design2.sim.list
}

n_grid <- seq(100,500,50)

design.2.sim.by.n.delta.2 <- vector("list", length = length(n_grid))

set.seed(87643)
for (i in 1:length(n_grid)){
  design.2.sim.by.n.delta.2[[i]] <- Design2Simulate(n_grid[i], delta, n_sim)
}
