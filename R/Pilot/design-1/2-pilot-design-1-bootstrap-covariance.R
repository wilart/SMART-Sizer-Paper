#Compute the covariance matrix
Design1ComputeCovAIPW <- function(n, delta, dat){
  #Computes AIPW Covariance for a given design 1 SMART data set
  #Source: ADHD simulation code from Prof. Ashkan Ertefaie
  
  dat$ID<-1:n
  dat.R1<-dat.R2<-dat[dat$r==1,];
  dat.R1$a2<-1; dat.R2$a2<- -1;
  rep.dat<-rbind(dat.R1, dat.R2, dat[dat$r==0,]); #Create rep data frame for estimation purposes
  w<-rep(1,dim(rep.dat)[1])
  w[rep.dat$s==1]<-2
  y<-rep.dat$y
  
  a1<-dat$a1
  a2<-dat$a2
  y2<-dat$y
  
  den<-num<-0
  
  #den2<-num2<-0; 
  meat<-meat2<-meat3<-meat4<-meat5<-meat6<-meat7<-meat8<-meat9<-meat10<-meat11<-meat12<-meat13<-bread2<-0
  
  fm<-glm(y~a1+a2, data=rep.dat,weight=w); summary(fm);
  
  cf<-coef(fm)
  
  AI.m<-matrix(c(1,1,-1,1,1,-1,-1,-1),ncol=2,byrow=T) #Contrast matrix (D)
  
  for(k in 1:4){
    a1i <- AI.m[k,1]
    a2j <- AI.m[k,2]
    ####### BELOW IS FOR POINT ESTIMATE
    cof <- coef(lm(y~o11+o12+o21+o22+a1+a1*o11+a2,data= rep.dat,weight=w) )
    
    des<-cbind(1,rep(a1i, n),rep(a2j,n))
    
    phi2<-cof[1]+cof[2]*dat$o11+cof[3]*dat$o12+ cof[4]*dat$o21+cof[5]*dat$o22+cof[6]*a1i+1*cof[7]*a2j+cof[8]*a1i*dat$o11
    
    phi22<-cof[1]+cof[2]*rep.dat$o11+cof[3]*rep.dat$o12+ cof[4]*rep.dat$o21+cof[5]*rep.dat$o22+cof[6]*rep.dat$a1+1*cof[7]*a2j+cof[8]*rep.dat$a1*rep.dat$o11
    
    cof2<-coef(lm(phi22~o11+o12+a1+a1*o11,data= rep.dat,weight=w) )
    phi1<-cof2[1]+cof2[2]*dat$o11+cof2[3]*dat$o12+cof2[4]*a1i+cof2[5]*a1i*dat$o11
    
    deni<-t(des)%*%(des)
    numi<- (1)*t(dat$r*des)%*%((y2-phi2+1*cof[7]*a2j)*((a1==a1i)/mean(a1==a1i)))+t(1*(1-dat$r)*des)%*%((y2-phi2)*((a1==a1i&a2==a2j)/mean(a1==a1i&a2==a2j)))+t(dat$r*des)%*%((phi2-phi1-1*cof[7]*a2j)*((a1==a1i)/mean(a1==a1i)))+t((1-dat$r)*des)%*%((phi2-phi1)*((a1==a1i)/mean(a1==a1i)))+t(1*des)%*%(phi1)
    
    den<-den+deni
    num<-num+numi
    #########################END of POINT ESTIMATION
    #######################VARIANCE ESTIMATION
    meat2 <- meat2+t(dat$r*des)%*%diag((y2-phi2+1*cof[7]*a2j)*((a1==a1i)/mean(a1==a1i)))
    meat3 <- meat3+t(1*(1-dat$r)*des)%*%diag((y2-phi2)*((a1==a1i&a2==a2j)/mean(a1==a1i&a2==a2j)))
    meat4 <- meat4+t(dat$r*des)%*%diag((phi2-phi1-1*cof[7]*a2j)*((a1==a1i)/mean(a1==a1i)))
    meat5 <- meat5+t((1-dat$r)*des)%*%diag((phi2-phi1)*((a1==a1i)/mean(a1==a1i)))
    meat6 <- meat6+t(des)%*%diag(as.vector(phi1-des%*% cf))
    bread2 <- bread2+t(des)%*%(des)
    
  }
  meat7 <- meat2+meat3+meat4+meat5+meat6		
  M <- meat7%*%t(meat7)
  
  DR.cov <- solve(bread2)%*%M%*%solve(bread2)*n
  ThetaCov <- cbind(rep(1,4),AI.m) %*% (DR.cov %*% t(cbind(rep(1,4),AI.m)))
  ThetaCov
}




set.seed(13954)
library(doParallel)
no_cores <- detectCores()
cl <- makeCluster(no_cores)
registerDoParallel(cl)

design.1.pilot.cov.bootstrap.list <- foreach(j = 1:1000) %dopar% {
    Design1ComputeCovAIPW(n = 50, delta = 0.25, design.1.pilot.bootstrap[[j]])
  }

stopCluster(cl)


library(Matrix)
load("~/../Dropbox/SMART-R/Rda/design.1.cov.true.list.rda")

design.1.cov.true <- nearPD(apply(simplify2array(design.1.cov.true.list),1:2, mean))$mat

design.1.pilot.cov.bootstrap.known.var <- vector("list", length = 1000)
for (j in 1:1000) {
  design.1.pilot.cov.bootstrap.known.var[[j]] <- sqrt(diag(diag(design.1.cov.true),4)) %*% cov2cor(design.1.pilot.cov.bootstrap.list[[j]]) %*% sqrt(diag(diag(design.1.cov.true),4))
}




design.1.pilot.sample.cov <- Design1ComputeCovAIPW(n = 50, delta = 0.25, design.1.pilot[[1]])
design.1.pilot.sample.known.var.cov <- sqrt(diag(diag(design.1.cov.true),4)) %*% cov2cor(design.1.pilot.sample.cov) %*% sqrt(diag(diag(design.1.cov.true),4))
