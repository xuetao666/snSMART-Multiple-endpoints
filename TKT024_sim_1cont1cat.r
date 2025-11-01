rm(list = ls())
library(R2WinBUGS)
library(rjags)
library(R2jags)
library(dplyr)
library(Matrix)


set.seed(19)


######################################################################################################################  
######################################################################################################################
#Generating data    
#####################################################################################################################



trial.sim.ran <- function(mu.true.stage1      = NULL, 
                          mu.true.stage2      = NULL, 
                          marginal.cat.p      = NULL,
                          n                   = NULL, 
                          cov.mat             = NULL,
                          sd.outcome2         = NULL) 
{
  # Total number of subjects
  n.length <- sum(n)
  max.trt.code <- length(n)
  
  # Stage I treatment assignment
  trt1 <- sort(rep(1:max.trt.code, times=n))
  
  # Stage II treatment assignment: randomly assign to low or high (2 or 3)
  
  trt2 <- sample(2:max.trt.code, n.length, replace = TRUE, prob=c(0.5,0.5))
  
  # Combine stage I and stage II treatments
  data.sim.1 <- data.frame(cbind(trt1, trt2))
  colnames(data.sim.1) <- c("trt1", "trt2")
  
  y <- matrix(0, nrow = n.length, ncol = 4)
  cov.mat.all <- bdiag(cov.mat, cov.mat)
  
  for (i in 1:n.length) {
    # Stage I and II outcome
    mu.y <- c(mu.true.stage1[1, trt1[i]],
              mu.true.stage1[2, trt1[i]],
              mu.true.stage2[1, trt2[i]],
              mu.true.stage2[2, trt2[i]]
    )
    
    y[i,1:4] <-MASS::mvrnorm(n = 1, mu = mu.y, Sigma = cov.mat.all)
    
  }
  
  data.sim.1 <- data.frame(cbind(data.sim.1,y))
  
  colnames(data.sim.1)<- c("trt1", "trt2", "V1.s1.cont", "V2.s1.cont","V1.s2.cont", "V2.s2.cont")
  data.sim.1$ID <- 1:nrow(data.sim.1)
  
  data.sim.s1 <- data.sim.1[,c("ID", "trt1", "V1.s1.cont", "V2.s1.cont")]
  data.sim.s2 <- data.sim.1[,c("ID", "trt2", "V1.s2.cont", "V2.s2.cont")]  
  
  
  #Defining cutpoints categorical variable 
  
  cut1=qnorm(marginal.cat.p[1],mean=mu.true.stage1[2,1],sd=sd.outcome2)
  cut2=qnorm(marginal.cat.p[2],mean=mu.true.stage1[2,1],sd=sd.outcome2)
  
  
  V2.s1.ord=ifelse(data.sim.1$V2.s1.cont<cut1,1,
                   ifelse(data.sim.1$V2.s1.cont>=cut1 & data.sim.1$V2.s1.cont<cut2,2,
                          ifelse(data.sim.1$V2.s1.cont>=cut2,3,NA)))
  
  V2.s2.ord=ifelse(data.sim.1$V2.s2.cont<cut1,1,
                   ifelse(data.sim.1$V2.s2.cont>=cut1 & data.sim.1$V2.s2.cont<cut2,2,
                          ifelse(data.sim.1$V2.s2.cont>=cut2,3,NA)))
  data.sim.f=cbind(data.sim.1,V2.s1.ord,V2.s2.ord)
  
  #Marginal.cat.stage1
  marginal.cat.stage1=marginal.cat.p
  for(j in 2:ncol(mu.true.stage1)){
    temp1=pnorm(cut1,mu.true.stage1[2,j],sd.outcome2)
    temp2=pnorm(cut2,mu.true.stage1[2,j],sd.outcome2)
    marginal.cat.stage1=rbind(marginal.cat.stage1,c(temp1,temp2))
  }
  
  return(list(data.sim.f=data.sim.f,marginal.cat.stage1=marginal.cat.stage1))
  
}

#Generating data under exchangeable and non-exchangeable scenarios

# True standard deviation values

sd.outcome1 <- 5
sd.outcome2 <- 2
rho.outcome12 <- 0.9

cov.matrix.true <- rbind(c(sd.outcome1^2, rho.outcome12*sd.outcome1*sd.outcome2),
                         c(rho.outcome12*sd.outcome1*sd.outcome2, sd.outcome2^2))



##True parameters

sample.n<-30
# Define sample sizes for each arm
n <- c(sample.n, sample.n, sample.n)  # Number of participants in each treatment arm

#Exchangeable scenario

# True mean values for each treatment in Stage I
mu.true.stage1 <- rbind(c(1.9, 1.9, 1.9),  
                        c(-3, -3, -3))    

mu.true.stage2 <-  rbind(c(1.9, 1.9, 1.9),  
                         c(-3, -3, -3)) 

marginal.cat.p<-c(0.33,0.66)


#Simulate data:


sim.result<- trial.sim.ran(mu.true.stage1= mu.true.stage1, 
                          mu.true.stage2= mu.true.stage2,
                          marginal.cat.p=marginal.cat.p,
                          n             = n, 
                          cov.mat       = cov.matrix.true,
                          sd.outcome2         = sd.outcome2) 
data.sim=sim.result$data.sim.f
marginal.cat.s1=sim.result$marginal.cat.stage1


######################################################################################################################

#Codes for models: bivariate model (parameter and likelihood) for BHM, bivariate model for BJSM 

#####################################################################################################################

BHM_1cont1cat.bugs<- function(){
  
  # Exchangeability Priors
  
  for (k in 1:3) {
    for (j in 1:2) {
      # between stage heterogeneity prior follows half normal distribution
      prec[j,k] <- pow(m.sd.prior[j], -2)
      tau[j,k] ~ dnorm(0, prec[j,k])%_% I(0,)
      
      # Exchangeable mean: gamma follows normal (mean, sd square)
      
      tau.mu.ex[j,k] <- pow(Prior.mu.ex[j,k,2], -2)
      gamma[j,k] ~ dnorm(Prior.mu.ex[j,k,1], tau.mu.ex[j,k])
      
    }
    
    # Latent variable for exchangeability 
    Z.exch[k,1] ~ dbin(p.exch[k, 1], 1)
    
    # Priors for rho
    rho.mu[k,1] ~ dunif(Prior.rho.mu[1], Prior.rho.mu[2])
  }
  
  
  for(k in 1:3){
    
    tau.mat.inv[k,1,1]<-pow(tau[1,k],2)
    tau.mat.inv[k,2,2]<-pow(tau[2,k],2)
    tau.mat.inv[k,1,2]<-rho.mu[k,1]*tau[1,k]*tau[2,k]
    tau.mat.inv[k,2,1]<-tau.mat.inv[k,1,2]
    tau.mat[k,1:2,1:2] <- inverse(tau.mat.inv[k,,])
    
    mu.ex[1:2,k] ~ dmnorm(gamma[1:2,k], tau.mat[k,1:2,1:2])
  }
  #NEX priors  
  
  for(k in 1:3){    
    for(s in 1:2){  
      v.mat.inv[k,s,1,1]<- pow(Prior.v.nex[k,s,1], 2)
      v.mat.inv[k,s,2,2]<- pow(Prior.v.nex[k,s,2], 2)
      v.mat.inv[k,s,1,2]<- rho.mu[k,1]*sqrt(v.mat.inv[k,s,1,1])*sqrt(v.mat.inv[k,s,2,2])
      v.mat.inv[k,s,2,1]<- v.mat.inv[k,s,1,2]
      v.mat[k,s,1:2,1:2] <- inverse(v.mat.inv[k,s,,])
      
      mu.nex[1:2,s,k] ~ dmnorm(Prior.mu.nex[k,s,1:2], v.mat[k,s,1:2,1:2])
    }
  }
  
  #Latent mixture assignment for ex-nex
  
  for(s in 1:2){
    for(k in 2:3){
      mu[1,s,k]<-Z.exch[k,1]*mu.ex[1,k]+(1-Z.exch[k,1])*mu.nex[1,s,k]
      mu[2,s,k]<-Z.exch[k,1]*mu.ex[2,k]+(1-Z.exch[k,1])*mu.nex[2,s,k]
    }
  }
  mu[1,1,1]<-mu.nex[1,1,1]
  mu[2,1,1]<-mu.nex[2,1,1]
  
  for(k in 1:3){ 
    sigma.1[k] ~ dunif(Prior.sigma.1[1], Prior.sigma.1[2])
    sigma.2[k] ~ dunif(Prior.sigma.2[1], Prior.sigma.2[2])
    prec.sigma1[k]<-pow(sigma.1[k],-2)
  }
  
  ##############################################################################
  #Likelihood variance-covariance 
  ##############################################################################
  
  for(k in 1:3){
    sigma.con[k]<-pow(sigma.2[k],2)*(1-pow(rho.y[k,1],2))
    prec.sigma.con[k]<-pow(sigma.con[k],-1)
  }
  
  
  for(k in 1:3){
    
    rho.y[k,1]~ dunif(Prior.rho.y[1], Prior.rho.y[2])
  }
  
  
  
  ##############################################################################
  #Likelihood calculation
  ##############################################################################
  for (i in 1:Nobs) {
    
    #Stage I
    mu.y1[i,1] <- mu[1,1,trt1[i]]  #outcome 1
    yy[i,1] ~ dnorm(mu.y1[i,1],prec.sigma1[trt1[i]])
    mu.z[i,1] <-mu[2,1,trt1[i]]+rho.y[trt1[i],1]*sigma.2[trt1[i]]*(yy[i,1]-mu.y1[i,1])/sigma.1[trt1[i]]
    z[i,1] ~ dnorm(mu.z[i,1],prec.sigma.con[trt1[i]])
    
    
    
    #Stage II
    mu.y2[i,1]<- mu[1,2,trt2[i]] #outcome 1
    yy[i,2] ~ dnorm(mu.y2[i,1],prec.sigma1[trt2[i]])
    mu.z[i,2] <-mu[2,2,trt2[i]]+rho.y[trt2[i],1]*sigma.2[trt2[i]]*(yy[i,2]-mu.y2[i,1])/sigma.1[trt2[i]]
    z[i,2] ~ dnorm(mu.z[i,2],prec.sigma.con[trt2[i]])
    
    
    #Transfer latented variable for outcome 2: Allows 
    
    y.cat[i,1]~dcat(pr.1[i,1:totcat1]) #Stage 1
    y.cat[i,2]~dcat(pr.2[i,1:totcat2]) #Stage 2
    
    for (r1 in 1:(totcat1-1)) {
      s.1[i,r1]<-phi((lambda.1[trt1[i],r1]- z[i,1])/sigma.2[trt1[i]])
    }
    
    for (rr in 1:(totcat2-1)) {
      s.2[i,rr]<-phi((lambda.2[trt2[i], rr]- z[i,2])/sigma.2[trt2[i]])
    }
    
    pr.1[i,1]<-s.1[i,1]
    pr.2[i,1]<-s.2[i,1]
    
    for(l in 2:(totcat1-1)) {
      pr.1[i,l]<-max(s.1[i,l]-s.1[i,l-1], 0.0001)
    }
    
    for(ll in 2:(totcat2-1)) {
      pr.2[i,ll]<-max(s.2[i,ll]-s.2[i,ll-1], 0.0001)
    }
    
    pr.1[i,totcat1]<- (1-s.1[i,totcat1-1])
    pr.2[i,totcat2]<- (1-s.2[i,totcat2-1])   
    
  }
  
  # lambda.11[1]<-0	
  for(k in 1:max.trt1){
    lambda.11[k, 1]<-0	
    for(ii in 2:(totcat1-1)){	
      lambda.11[k, ii]~ dnorm(0,1)%_% I(0,)	
    }
    
    lambda.1[k,1:2] <- sort(lambda.11[k,1:2])
  }
  
  for(kk in 1:max.trt2){
    
    lambda.22[kk, 1]<- 0	
    for(ii in 2:(totcat2-1)){	
      lambda.22[kk, ii]~ dnorm(0,1)%_% I(0,)	
    }
    
    lambda.2[kk,1:2] <- sort(lambda.22[kk,1:2])
  }
  
  #Adding global inference
  mu.diff[1,2] <-mu[1,1,2]-mu[1,1,1]
  mu.diff[1,3] <-mu[1,1,3]-mu[1,1,1]
  
}    


#Model

model.BHM <- BHM_1cont1cat.bugs

m.sd.prior<-c(2,1)

# m.sd.prior <- c(0.25, 0.25, 0.25)
#Update prior.mu.ex
Prior.mu.ex <- array(sort(rep(c(0,10), times=6)), dim=c(2,3,2)) 
# Prior.mu.ex[, , 1]<-mean.prior
# Prior.mu.ex[, , 2]<-sd.prior

p.exch <- cbind(c(0,0.5,0.5), c(0,0.5,0.5))
Prior.rho.mu <- c(-1,1)
Prior.mu.nex <- array(rep(c(0,0), times=6), dim=c(3,2,2)) 
# Prior.mu.nex[, , 1]<-mean.prior
# Prior.mu.nex[, , 2]<-sd.prior

Prior.v.nex <- array(sort(rep(c(5,1), times=6),decreasing =TRUE), dim=c(3,2,2))
Prior.v.nex[1,1:2,1:2]<-sqrt(10)

Prior.sigma.1 <- c(0,10)
Prior.sigma.2 <- c(0,10)
Prior.rho.y <- c(-1,1)


MCMC  <- c(40000, 50000, 2, 2)


Nobs <- sum(n)
trt1 <- data.sim$trt1
trt2 <- data.sim$trt2



yy=cbind(data.sim$V1.s1.cont,data.sim$V1.s2.cont)

y.cat=cbind(data.sim$V2.s1.ord,data.sim$V2.s2.ord)
totcat1<- 3
totcat2<- 3
max.trt1=3
max.trt2=3
#Model 1: BHM with correlation 

data.BHM <- list("Nobs", "trt1", "trt2", "yy","y.cat","totcat1","totcat2",
                 "m.sd.prior", "Prior.mu.ex", "p.exch",
                 "Prior.rho.mu", "Prior.mu.nex", "Prior.v.nex", 
                 "Prior.sigma.1", "Prior.sigma.2", "Prior.rho.y",
                 "max.trt1","max.trt2")

initsfun.BHM = function(i){
  list(
    tau = matrix(rexp(6,10), nrow=2, ncol=3),
    gamma = matrix(rnorm(6,0,10), nrow=2, ncol=3),
    rho.mu = cbind(NULL, runif(3, -1, 1)),
    Z.exch = cbind(NULL, rep(0,3)),
    mu.ex = matrix(rnorm(6, 0, 10), nrow=2, ncol=3),
    mu.nex= array(rnorm(12,0,10),dim= c(2, 2, 3)),
    sigma.1 = rexp(3, 1),
    sigma.2 = rexp(3, 1),
    rho.y = cbind(NULL,runif(3, -1, 1))
  )
}

inits.BHM <- lapply(rep(1, MCMC[3]),  initsfun.BHM)

param.BHM <- c("mu", "mu.diff","rho.y","sigma.1","sigma.2","rho.mu","pr.1","pr.2","lambda.1","lambda.2")

BHM.fit<- jags(
  data = data.BHM,
  model = model.BHM,
  parameters.to.save = param.BHM,
  n.burnin = MCMC[1],
  n.iter = MCMC[2],
  n.chains = MCMC[3],
  n.thin = MCMC[4],
  inits = inits.BHM,
  DIC = TRUE,
  progress.bar= "none",
  jags.seed = array.id 
)

# ##########################################################################################  
# #########################################################################################  
# ##(BJSM with  bivariate model for parameter and likelihood)
# ##########################################################################################
# 
BJSM_1cont1cat.bugs<- function(){
  # Priors definition
  
  
  for (k in 1:3) {
    
    tau.mat.inv[k,1,1]<-pow(Prior.mu.ex[1,k,2],2)
    tau.mat.inv[k,2,2]<-pow(Prior.mu.ex[2,k,2],2)
    tau.mat.inv[k,1,2]<-rho.mu[k,1]*Prior.mu.ex[1,k,2]*Prior.mu.ex[2,k,2]
    tau.mat.inv[k,2,1]<-tau.mat.inv[k,1,2]
    tau.mat[k,1:2,1:2] <- inverse(tau.mat.inv[k,,])
    
    
    mu[1:2, k] ~ dmnorm(Prior.mu.ex[1:2,k,1], tau.mat[k,1:2,1:2])
    
  }
  
  
  # Priors for rho
  ##############################################################################
  for (k in 1:3){
    rho.mu[k,1] ~ dunif(-1, 1)
  }
  
  
  # Priors for alpha
  ##############################################################################
  for (j in 1:2){
    
    alpha[j] ~ dnorm(0,0.1)
  }
  
  # Joint model for mu
  for(j in 1:2){
    for(k in 1:3){
      mu.vec[j,1,k]<-mu[j,k]
      mu.vec[j,2,k]<-mu[j,k]+alpha[j]
    }
  }
  
  for(k in 1:3){ 
    sigma.1[k] ~ dunif(Prior.sigma.1[1], Prior.sigma.1[2])
    sigma.2[k] ~ dunif(Prior.sigma.2[1], Prior.sigma.2[2])
    prec.sigma1[k]<-pow(sigma.1[k],-2)
  }
  
  ##############################################################################
  #Likelihood variance-covariance 
  ##############################################################################
  
  for(k in 1:3){
    sigma.con[k]<-pow(sigma.2[k],2)*(1-pow(rho.y[k,1],2))
    prec.sigma.con[k]<-pow(sigma.con[k],-1)
  }
  
  
  for(k in 1:3){
    
    rho.y[k,1]~ dunif(Prior.rho.y[1], Prior.rho.y[2])
  }
  
  ##############################################################################
  #Likelihood calculation
  ##############################################################################
  for (i in 1:Nobs) {
    
    #Stage I
    mu.y1[i,1] <- mu.vec[1,1,trt1[i]]  #outcome 1
    yy[i,1] ~ dnorm(mu.y1[i,1],prec.sigma1[trt1[i]])
    mu.z[i,1] <-mu.vec[2,1,trt1[i]]+rho.y[trt1[i],1]*sigma.2[trt1[i]]*(yy[i,1]-mu.y1[i,1])/sigma.1[trt1[i]]
    z[i,1] ~ dnorm(mu.z[i,1],prec.sigma.con[trt1[i]])
    
    
    
    #Stage II
    mu.y2[i,1]<- mu.vec[1,2,trt2[i]] #outcome 1
    yy[i,2] ~ dnorm(mu.y2[i,1],prec.sigma1[trt2[i]])
    mu.z[i,2] <-mu.vec[2,2,trt2[i]]+rho.y[trt2[i],1]*sigma.2[trt2[i]]*(yy[i,2]-mu.y2[i,1])/sigma.1[trt2[i]]
    z[i,2] ~ dnorm(mu.z[i,2],prec.sigma.con[trt2[i]])
    
    
    #Transfer latented variable for outcome 2: Allows 
    
    y.cat[i,1]~dcat(pr.1[i,1:totcat1]) #Stage 1
    y.cat[i,2]~dcat(pr.2[i,1:totcat2]) #Stage 2
    
    for (r1 in 1:(totcat1-1)) {
      s.1[i,r1]<-phi((lambda.1[trt1[i],r1]- z[i,1])/sigma.2[trt1[i]])
    }
    
    for (rr in 1:(totcat2-1)) {
      s.2[i,rr]<-phi((lambda.2[trt2[i], rr]- z[i,2])/sigma.2[trt2[i]])
    }
    
    pr.1[i,1]<-s.1[i,1]
    pr.2[i,1]<-s.2[i,1]
    
    for(l in 2:(totcat1-1)) {
      pr.1[i,l]<-max(s.1[i,l]-s.1[i,l-1], 0.0001)
    }
    
    for(ll in 2:(totcat2-1)) {
      pr.2[i,ll]<-max(s.2[i,ll]-s.2[i,ll-1], 0.0001)
    }
    
    pr.1[i,totcat1]<- (1-s.1[i,totcat1-1])
    pr.2[i,totcat2]<- (1-s.2[i,totcat2-1])   
    
  }
  
  # lambda.11[1]<-0	
  for(k in 1:max.trt1){
    lambda.11[k, 1]<-0	
    for(ii in 2:(totcat1-1)){	
      lambda.11[k, ii]~ dnorm(0,1)%_% I(0,)	
    }
    
    lambda.1[k,1:2] <- sort(lambda.11[k,1:2])
  }
  
  for(kk in 1:max.trt2){
    
    lambda.22[kk, 1]<- 0	
    for(ii in 2:(totcat2-1)){	
      lambda.22[kk, ii]~ dnorm(0,1)%_% I(0,)	
    }
    
    lambda.2[kk,1:2] <- sort(lambda.22[kk,1:2])
  }
  
  #Adding global inference
  mu.diff[1,2] <-mu.vec[1,1,2]-mu.vec[1,1,1]
  mu.diff[1,3] <-mu.vec[1,1,3]-mu.vec[1,1,1]
  
}




#Model

model.BJSM <- BJSM_1cont1cat.bugs
data.BJSM <-  list("Nobs", "trt1", "trt2", "yy","y.cat","totcat1","totcat2","Prior.rho.y",
                   "Prior.sigma.1", "Prior.sigma.2","Prior.mu.ex",
                   "max.trt1","max.trt2")

initsfun.BJSM = function(i){
  list(
    mu = matrix(rnorm(6, 0,10), nrow=2, ncol=3),
    rho.mu = cbind(NULL,runif(3, -1, 1)),
    alpha = rnorm(2, 0,10),
    sigma.1 = rexp(3, 1),
    sigma.2 = rexp(3, 1),
    rho.y = cbind(NULL,runif(3, -1, 1))
  )
}

inits.BJSM <- lapply(rep(1, MCMC[3]),  initsfun.BJSM)

param.BJSM <- c("mu", "mu.diff","rho.y","rho.mu","sigma.1","sigma.2","alpha","mu.vec","pr.1","pr.2","lambda.1")

BJSM.fit<- jags(
  data = data.BJSM,
  model = model.BJSM,
  parameters.to.save = param.BJSM,
  n.burnin = MCMC[1],
  n.iter = MCMC[2],
  n.chains = MCMC[3],
  n.thin = MCMC[4],
  inits = inits.BJSM,
  DIC = TRUE,
  progress.bar= "none",
  jags.seed = array.id
)




#Summary part:
#-------------------------------------------------------------------------------
#BHM
#-------------------------------------------------------------------------------

true.cat=cbind(marginal.cat.s1[,1],marginal.cat.s1[,2]-marginal.cat.s1[,1],1-marginal.cat.s1[,2])

mean.true <- c(mu.true.stage1[1,2:3]-mu.true.stage1[1,1],
               mu.true.stage1[1,1],true.cat[1,1]-true.cat[2,1],true.cat[1,1]-true.cat[3,1],true.cat[1,],true.cat[2,],true.cat[3,])


get_probs_one_treatment <- function(data, j, eps = 1e-12) {
  mu  <- data[, sprintf("mu[2,1,%d]", j)]
  l1  <- data[, sprintf("lambda.1[%d,1]", j)]  
  l2  <- data[, sprintf("lambda.1[%d,2]", j)]  
  sd  <- data[, sprintf("sigma.2[%d]",j)]
  
  s1 <- pnorm(l1 - mu,sd=sd)
  s2 <- pnorm(l2 - mu,sd=sd)
  
  p1 <- pmax(s1, eps)
  p2 <- pmax(s2 - s1, eps)
  p3 <- pmax(1 - s2, eps)
  
  den <- p1 + p2 + p3
  cbind(cat1 = p1/den, cat2 = p2/den, cat3 = p3/den)
}

{
  
  #Output
  #--------------------------------------------------------------------------------------------------------------
  sum=round(BHM.fit$BUGSoutput$summary,2)
  BHM=BHM.fit$BUGSoutput$sims.matrix
  # Compute mean probabilities for each category at Stage 1
  
  
  
  pr.1_samples <- BHM.fit$BUGSoutput$sims.list$pr.1
  
  
  pr.1=NULL
  id=1:(sample.n*3)
  for(trti in 1:3){
    i_vec <- id[data.sim$trt1 == trti]
    pr.1<-cbind(pr.1,apply(pr.1_samples[,data.sim$trt1 == trti,],c(1,3),mean))
  }
  
  #USE lambda:
  
  P_j1 <- get_probs_one_treatment(BHM, 1)
  P_j2 <- get_probs_one_treatment(BHM, 2)
  P_j3 <- get_probs_one_treatment(BHM, 3)
  
  lambda.p=cbind(P_j1,P_j2,P_j3)
  colMeans(lambda.p)
  colMeans(pr.1)
  
  
  
  BHM=cbind(BHM[,"mu.diff[1,2]"],BHM[,"mu.diff[1,3]"],
            BHM[,"mu[1,1,1]"],lambda.p[,1]-lambda.p[,4],lambda.p[,1]-lambda.p[,7],lambda.p)
  BHM.mean=colMeans(BHM)
  #Bias
  BHM.Bias <- round(BHM.mean-mean.true,2)
  
  #Coverage:
  BHM.lower <- apply(BHM, 2, quantile, probs = 0.025)
  BHM.upper <- apply(BHM, 2, quantile, probs = 0.975)
  
  # Check if true value is inside the interval
  BHM.coverage <- (mean.true >= BHM.lower) & (mean.true <= BHM.upper)
  
  
  #Global inference:
  BHM.mu.diff=BHM.fit$BUGSoutput$sims.matrix[,c("mu.diff[1,2]",
                                                "mu.diff[1,3]")]
  
  ###Either low or high doses have significant different than placebo in either stages:
  ##Outcome 1
  BHM.out1.inf=(sum(BHM.mu.diff[,"mu.diff[1,2]"]>0 | BHM.mu.diff[,"mu.diff[1,3]"]>0)/nrow(BHM.mu.diff))>0.99375
  BHM.out1.l.inf=(sum(BHM.mu.diff[,"mu.diff[1,2]"]>0)/nrow(BHM.mu.diff))>0.9875
  BHM.out1.h.inf=(sum(BHM.mu.diff[,"mu.diff[1,3]"]>0)/nrow(BHM.mu.diff))>0.9875
  ##Outcome 2
  #low
  BHM.out2.l.inf=sum(lambda.p[,1]>lambda.p[,4])/nrow(lambda.p)>0.9875
  BHM.out2.h.inf=sum(lambda.p[,1]>lambda.p[,7])/nrow(lambda.p)>0.9875
  BHM.out2.inf=sum(lambda.p[,1]>lambda.p[,4] | lambda.p[,1]>lambda.p[,7])/nrow(lambda.p)>0.99375
  #Either
  BHM.either.l.inf=sum(BHM.mu.diff[,"mu.diff[1,2]"]>0 | lambda.p[,1]>lambda.p[,4])/nrow(lambda.p)>0.99375
  BHM.either.h.inf=sum(BHM.mu.diff[,"mu.diff[1,3]"]>0 | lambda.p[,1]>lambda.p[,7])/nrow(lambda.p)>0.99375
  BHM.either.inf=sum(BHM.mu.diff[,"mu.diff[1,2]"]>0 | BHM.mu.diff[,"mu.diff[1,3]"]>0 |
                       lambda.p[,1]>lambda.p[,4] | lambda.p[,1]>lambda.p[,7])/nrow(lambda.p)>0.99844
  
  
  
  BHM.inf=c(BHM.out1.l.inf,BHM.out1.h.inf,BHM.out1.inf,BHM.out2.l.inf,BHM.out2.h.inf,BHM.out2.inf,
            BHM.either.l.inf,BHM.either.h.inf,BHM.either.inf)
  BHM.absbias=abs(BHM.Bias)
}


#-------------------------------------------------------------------------------
#BJSM
#-------------------------------------------------------------------------------
get_probs_one_treatment_bjsm <- function(data, j, eps = 1e-12) {
  mu  <- data[, sprintf("mu.vec[2,1,%d]", j)]
  l1  <- data[, sprintf("lambda.1[%d,1]", j)]  
  l2  <- data[, sprintf("lambda.1[%d,2]", j)]
  sd  <- data[, sprintf("sigma.2[%d]", j)]
  
  s1 <- pnorm(l1 - mu,sd=sd)
  s2 <- pnorm(l2 - mu,sd=sd)
  
  p1 <- pmax(s1, eps)
  p2 <- pmax(s2 - s1, eps)
  p3 <- pmax(1 - s2, eps)
  
  den <- p1 + p2 + p3
  cbind(cat1 = p1/den, cat2 = p2/den, cat3 = p3/den)
}
{
  
  #Output
  #--------------------------------------------------------------------------------------------------------------
  sum=round(BJSM.fit$BUGSoutput$summary,2)
  BJSM=BJSM.fit$BUGSoutput$sims.matrix
  # Compute mean probabilities for each category at Stage 1
  pr.1_samples <- BJSM.fit$BUGSoutput$sims.list$pr.1
  
  
  pr.1=NULL
  id=1:(sample.n*3)
  for(trti in 1:3){
    i_vec <- id[data.sim$trt1 == trti]
    pr.1<-cbind(pr.1,apply(pr.1_samples[,data.sim$trt1 == trti,],c(1,3),mean))
  }
  
  
  P_j1 <- get_probs_one_treatment_bjsm(BJSM, 1)
  P_j2 <- get_probs_one_treatment_bjsm(BJSM, 2)
  P_j3 <- get_probs_one_treatment_bjsm(BJSM, 3)
  
  lambda.p=cbind(P_j1,P_j2,P_j3)
  colMeans(pr.1)
  colMeans(lambda.p)
  
  
  
  
  BJSM=cbind(BJSM[,"mu.diff[1,2]"],BJSM[,"mu.diff[1,3]"],
             BJSM[,"mu[1,1]"],lambda.p[,1]-lambda.p[,4],lambda.p[,1]-lambda.p[,7],lambda.p)
  BJSM.mean=colMeans(BJSM)
  #Bias
  BJSM.Bias <- round(BJSM.mean-mean.true,2)
  
  #Coverage:
  BJSM.lower <- apply(BJSM, 2, quantile, probs = 0.025)
  BJSM.upper <- apply(BJSM, 2, quantile, probs = 0.975)
  
  # Check if true value is inside the interval
  BJSM.coverage <- (mean.true >= BJSM.lower) & (mean.true <= BJSM.upper)
  
  
  #Global inference:
  BJSM.mu.diff=BJSM.fit$BUGSoutput$sims.matrix[,c("mu.diff[1,2]",
                                                  "mu.diff[1,3]")]
  
  ###Either low or high doses have significant different than placebo in either stages:
  ##Outcome 1
  BJSM.out1.inf=(sum(BJSM.mu.diff[,"mu.diff[1,2]"]>0 | BJSM.mu.diff[,"mu.diff[1,3]"]>0)/nrow(BJSM.mu.diff))>0.99375
  BJSM.out1.l.inf=(sum(BJSM.mu.diff[,"mu.diff[1,2]"]>0)/nrow(BJSM.mu.diff))>0.9875
  BJSM.out1.h.inf=(sum(BJSM.mu.diff[,"mu.diff[1,3]"]>0)/nrow(BJSM.mu.diff))>0.9875
  ##Outcome 2
  #low
  BJSM.out2.l.inf=sum(lambda.p[,1]>lambda.p[,4])/nrow(lambda.p)>0.9875
  BJSM.out2.h.inf=sum(lambda.p[,1]>lambda.p[,7])/nrow(lambda.p)>0.9875
  BJSM.out2.inf=sum(lambda.p[,1]>lambda.p[,4] | lambda.p[,1]>lambda.p[,7])/nrow(lambda.p)>0.99375
  #Either
  BJSM.either.l.inf=sum(BJSM.mu.diff[,"mu.diff[1,2]"]>0 | lambda.p[,1]>lambda.p[,4])/nrow(lambda.p)>0.99375
  BJSM.either.h.inf=sum(BJSM.mu.diff[,"mu.diff[1,3]"]>0 | lambda.p[,1]>lambda.p[,7])/nrow(lambda.p)>0.99375
  BJSM.either.inf=sum(BJSM.mu.diff[,"mu.diff[1,2]"]>0 | BJSM.mu.diff[,"mu.diff[1,3]"]>0 |
                        lambda.p[,1]>lambda.p[,4] | lambda.p[,1]>lambda.p[,7])/nrow(lambda.p)>0.99844
  
  
  
  BJSM.inf=c(BJSM.out1.l.inf,BJSM.out1.h.inf,BJSM.out1.inf,BJSM.out2.l.inf,BJSM.out2.h.inf,BJSM.out2.inf,
             BJSM.either.l.inf,BJSM.either.h.inf,BJSM.either.inf)
  BJSM.absbias=abs(BJSM.Bias)
}
