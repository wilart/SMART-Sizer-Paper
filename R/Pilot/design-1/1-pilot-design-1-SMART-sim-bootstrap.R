#Simulate pilot SMART with 4 EDTRs and bootstrap
delta = 0.25

###
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

#c(32, 44, 200, 26) # based off 4 of 5 from https://methodology.psu.edu/ra/adap-inter/projects accessed 7-4-2018 https://www.ncbi.nlm.nih.gov/pubmed/25785788, https://clinicaltrials.gov/ct2/show/NCT01880814, https://clinicaltrials.gov/ct2/show/NCT02414074, https://www.tandfonline.com/doi/abs/10.1080/15374416.2015.1102069

set.seed(12376)
design.1.pilot <- Design1Simulate(n = 50, delta = 0.25, n_sim = 1)

design.1.pilot.bootstrap <- vector("list", length = 1)
  for (j in 1:1000) {
    design.1.pilot.bootstrap[[j]] <- design.1.pilot[[1]][sample(50,replace = TRUE),]
    
  }
