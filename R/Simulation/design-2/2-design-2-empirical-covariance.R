#Compute Empirical Covariance

Design2ComputeCovAIPW <- function(n, delta, dat) {
  emp.cov<-cov2<-0
  
  
  #for(i in 1:irep){
  probS2<-c(mean(dat$a1==-1&dat$a2==1),mean(dat$a1==-1&dat$a2==2),mean(dat$a1==-1&dat$a2==3),mean(dat$a1==-1&dat$a2==4))
  mean(dat$r)
  dat$ID<-1:n
  dat$a2[dat$a1==1]<-0
  dat$a2[dat$a1==-1&dat$r==1]<-4
  
  
  dat.R1<-dat.R2<-dat.R3<-dat.R4<-dat[dat$a1==-1&dat$r==1,];
  dat.R1$a2<-1; dat.R2$a2<- 2;
  dat.R3$a2<-3; dat.R4$a2<- 4;
  
  rep.dat<-rbind(dat.R1, dat.R2,dat.R3, dat.R4, dat[dat$a1==1|dat$r==0,]); dim(rep.dat);
  w<-rep(1,dim(rep.dat)[1])
  w[rep.dat$a1==-1&rep.dat$r==0]<-4
  
  
  y<-rep.dat$y
  z<-dat$a1==-1
  
  rep.dat$z<-rep.dat$a1==-1
  a2z<-cbind(rep.dat$z*I(rep.dat$a2==1),rep.dat$z*I(rep.dat$a2==2),rep.dat$z*I(rep.dat$a2==3),rep.dat$z*I(rep.dat$a2==2)*rep.dat$o21)
  cof<-coef(lm(y~o11+o12+o21+o22+I(a1==-1)+I(a1==-1)*o11+a2z,data= rep.dat,weight=w) )
  cof
  
  phi2<-cof[1]+cof[2]*dat$o11+cof[3]*dat$o12+ cof[4]*dat$o21+cof[5]*dat$o22+cof[6]*(dat$a1==-1)+1*cof[7]*z*I(dat$a2==1)+cof[8]*z*I(dat$a2==2)+cof[9]*z*I(dat$a2==3)+cof[10]*(dat$a1==-1)*dat$o11
  
  phi22<-cof[1]+cof[2]*rep.dat$o11+cof[3]*rep.dat$o12+ cof[4]*rep.dat$o21+cof[5]*rep.dat$o22+cof[6]*(rep.dat$a1==-1)+1*cof[7]*a2z[1]+cof[8]*a2z[2]+cof[9]*a2z[3]+cof[10]*(rep.dat$a1==-1)*rep.dat$o11
  
  
  cof2<-coef(lm(phi22~o11+o12+a1+a1*o11,data= rep.dat,weight=w) )
  phi1<-cof2[1]+cof2[2]*dat$o11+cof2[3]*dat$o12+cof2[4]*dat$a1+cof2[5]*dat$a1*dat$o11
  
  
  
  
  cbind(rep.dat$a2,rep.dat$r,rep.dat$a1)
  rep.dat$a2[(rep.dat$a1==1)]<-4
  
  
  a2z<-cbind(rep.dat$z*I(rep.dat$a2==1),rep.dat$z*I(rep.dat$a2==2),rep.dat$z*I(rep.dat$a2==3))
  
  fm<-glm(y~a1+a2z, data=rep.dat,weight=w); summary(fm);
  cf<-coef(fm)
  
  AI.m<-matrix(c(1,0,0,0,-1,0,0,0,-1,1,0,0,-1,0,1,0,-1,0,0,1),ncol=4,byrow=T)
  #AI.m<-matrix(c(0,0,0,0,1,0,0,0,1,0,1,0,1,0,0,1),ncol=4,byrow=T)
  AI.m2<-matrix(c(1,0,-1,4,-1,1,-1,2,-1,3),ncol=2,byrow=T)
  
  
  rdn<-dim(rep.dat)[1]
  y2<-dat$y
  
  a1<-dat$a1
  a2<-dat$a2
  
  den<-num<-0
  den2<-num2<-0; meat<-meat2<-meat3<-meat4<-meat5<-meat6<-meat7<-meat8<-meat9<-meat10<-meat11<-meat12<-meat13<-bread2<-0
  
  for(k in 2:5){
    
    a1i<-AI.m[k,1]
    a2j<-AI.m2[k,2]
    
    a21j<-AI.m[k,2]
    a22j<-AI.m[k,3]
    a23j<-AI.m[k,4]
    
    ####### BELOW IS FOR POINT ESTIMATE
    
    #    	 des<-cbind(1,rep(a1i, n),(a1i==-1)*rep(a21j,n),(a1i==-1)*rep(a22j,n),(a1i==-1)*rep(a23j,n));
    
    des<-cbind(1,rep(as.numeric(a1i==-1), n),(a1i==-1)*rep(a21j,n),(a1i==-1)*rep(a22j,n),(a1i==-1)*rep(a23j,n));
    #design matrix?
    
    phi2<-cof[1]+cof[2]*dat$o11+cof[3]*dat$o12+ cof[4]*dat$o21+cof[5]*dat$o22+cof[6]*(a1i==-1)+cof[7]*a21j+1*cof[8]*a22j+1*cof[9]*a23j+1*cof[10]*a22j*dat$o21+cof[11]*(a1i==-1)*dat$o11
    
    
    phi22<-cof[1]+cof[2]* rep.dat$o11+cof[3]*rep.dat$o12+ cof[4]*rep.dat$o21+cof[5]*rep.dat$o22+cof[6]*(rep.dat$a1==-1)+1*cof[7]*a21j+1*cof[8]*a22j+1*cof[9]*a23j+1*cof[10]*a22j*rep.dat$o21+cof[11]*(rep.dat$a1==-1)*rep.dat$o11
    
    
    
    
    cof2<-coef(lm(phi22~o11+o12+I(a1==-1)+I(a1==-1)*o11,data= rep.dat,weight=w) )
    phi1<-cof2[1]+cof2[2]*dat$o11+cof2[3]*dat$o12+cof2[4]*(a1i==-1)+cof2[5]*(a1i==-1)*dat$o11
    
    #####################################
    deni<-t(des)%*%(des)
    
    numi<-t(dat$r*des)%*%((y2-phi2)*(a1==-1)/mean(a1==-1))+t((1-dat$r)*des)%*%((y2-phi2)*((a1==-1&a2==a2j)/probS2[k-1]))+t(dat$r*des)%*%((phi2-phi1)*((a1==-1)/mean(a1==-1)))+t((1-dat$r)*des)%*%((phi2-phi1)*((a1==-1)/mean(a1==-1)))+t(des)%*%(phi1)
    
    den<-den+deni
    num<-num+numi
    
    #######################VARIANCE ESTIMATION
    meat2<-meat2+t(dat$r*des)%*%diag((y2-phi2)*(a1==-1)/mean(a1==-1))+t((1-dat$r)*des)%*%diag((y2-phi2)*((a1==-1&a2==a2j)/probS2[k-1]))+t(dat$r*des)%*%diag((phi2-phi1)*((a1==-1)/mean(a1==-1)))+t((1-dat$r)*des)%*%diag((phi2-phi1)*((a1==-1)/mean(a1==-1)))+t(des)%*%diag(as.vector(phi1-des%*% cf))
    
    bread2<-bread2+t(des)%*%(des)
    
    
  }
  
  
  des<-cbind(1,rep(0, n),0,0,0);
  den<-den+t(des)%*%((des))
  phi2<-cof[1]+cof[2]*dat$o11+cof[3]*dat$o12+ cof[4]*dat$o21+cof[5]*dat$o22+cof[6]*0+cof[11]*0*dat$o11
  
  #phi22<-cof[1]+cof[2]*rep.dat$o11+cof[3]*rep.dat$o12+ cof[4]*rep.dat$o21+cof[5]*rep.dat$o22+cof[6]*rep.dat$a1+1*cof[7]*a2z[1]+cof[8]*a2z[2]+cof[9]*a2z[3]+cof[11]*rep.dat$a1*rep.dat$o11
  phi22<-cof[1]+cof[2]*rep.dat$o11+cof[3]*rep.dat$o12+ cof[4]*rep.dat$o21+cof[5]*rep.dat$o22+cof[6]*(rep.dat$a1==-1)+cof[11]*(rep.dat$a1==-1)*rep.dat$o11
  
  
  cof2<-coef(lm(phi22~o11+o12+I(a1==-1)+I(a1==-1)*o11,data= rep.dat,weight=w) )
  phi1<-cof2[1]+cof2[2]*dat$o11+cof2[3]*dat$o12+cof2[4]*0+cof2[5]*0*dat$o11
  
  #cbind(dat$a1,dat$a2)
  
  num<-num+t(des)%*%((y2-phi2)*(a1==+1)/mean(a1==1))+t(des)%*%((phi2-phi1)*((a1==+1)/mean(a1==1)))+t(des)%*%(phi1)
  
  meat3<-t(des)%*%diag((y2-phi2)*(a1==+1)/mean(a1==1))+t(des)%*%diag((phi2-phi1)*((a1==+1)/mean(a1==1)))+t(des)%*%diag(as.vector(phi1-des%*% cf))
  
  meat7<-meat2+meat3
  M<-meat7%*%t(meat7)
  bread2<-bread2+t(des)%*%(des)
  
  
  emp.cov<-solve(bread2)%*%M%*%solve(bread2) # THIS ESTIMATES THE
  
  DR.cov<-emp.cov #AIPW
  c1<-c(1,0,0,0,0);c2<-c(1,1,0,0,0);c3<-c(1,1,1,0,0);c4<-c(1,1,0,1,0);c5<-c(1,1,0,0,1)
  cont<-matrix(c(c1,c2,c3,c4,c5),nrow=5,byrow=TRUE)
  V<-cont%*%DR.cov%*%t(cont)*n
  V
}

design.2.cov.by.n.delta.2 <- vector("list", length = length(n_grid))
library(doParallel)
no_cores <- detectCores()
cl <- makeCluster(no_cores)
registerDoParallel(cl)


design.2.cov.by.n.delta.2 <- foreach(i = 1:length(n_grid)) %:%
  foreach(j = 1:n_sim) %dopar% {
    Design2ComputeCovAIPW(n_grid[i], delta, design.2.sim.by.n.delta.2[[i]][[j]])
  }
stopCluster(cl)

#save(design.2.cov.by.n.delta.2, file = '~/../Dropbox/SMART-R/Rda/design.2.cov.by.n.delta.2.rda')
