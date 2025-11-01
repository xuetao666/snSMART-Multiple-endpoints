

library(R2WinBUGS)
library(rjags)
library(R2jags)
library(dplyr)
library(Matrix)


set.seed(1)


######################################################################################################################  
######################################################################################################################
#Generating data    
#####################################################################################################################


trial.sim.ran <- function(mu.true.stage1= NULL, 
                          mu.true.stage2= NULL, 
                          n             = NULL, 
                          cov.mat       = NULL) 
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
  
  data.sim.1 <- cbind(data.sim.1,y)
  return(data.sim.1)
}


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
mu.true.stage1.ex <- rbind(c(1.9, 1.9, 1.9),  
                           c(-3, -3, -3))    

mu.true.stage2.ex <- rbind(c(1.9, 1.9, 1.9),  
                           c(-3, -3, -3)) 
#Simulate data:


data.sim <- trial.sim.ran(mu.true.stage1= mu.true.stage1.ex, 
                          mu.true.stage2= mu.true.stage2.ex, 
                          n             = n, 
                          cov.mat       = cov.matrix.true) 


##########################################################################################  
##########################################################################################    
#(BHM with  bivariate model for both parameter and likelihood (M3))
##########################################################################################

BHM_2cont.bugs<- function(){
  
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
  }
  
  ##############################################################################
  #Likelihood variance-covariance 
  ##############################################################################
  for(s in 1:2){
    for(k in 1:3){
      
      sigma.lkd[s,k,1,1] <- pow(sigma.1[k], 2)
      sigma.lkd[s,k,2,2] <- pow(sigma.2[k], 2)
      sigma.lkd[s,k,1,2] <- rho.y[k,1]*sqrt(sigma.lkd[s,k,1,1])*sqrt(sigma.lkd[s,k,2,2])
      sigma.lkd[s,k,2,1] <- sigma.lkd[s,k,1,2]
      prec.lkd[s,k,1:2,1:2] <- inverse(sigma.lkd[s,k,,])
    }
  }
  
  for(k in 1:3){
    
    rho.y[k,1]~ dunif(Prior.rho.y[1], Prior.rho.y[2])
  }
  
  
  
  ##############################################################################
  #Likelihood calculation
  ##############################################################################
  for (i in 1:Nobs) {
    
    mu.y1[i,1] <- mu[1,1,trt1[i]]
    mu.y1[i,2] <- mu[2,1,trt1[i]]
    mu.y2[i,1] <- mu[1,2,trt2[i]]
    mu.y2[i,2] <- mu[2,2,trt2[i]]
    
    y.new[i,1:2] ~ dmnorm(mu.y1[i,1:2], prec.lkd[1,trt1[i],,])
    y.new[i,3:4] ~ dmnorm(mu.y2[i,1:2], prec.lkd[2,trt2[i],,])
  }
  
  #Adding global inference
  mu.diff[1,1,2] <-mu[1,1,2]-mu[1,1,1]
  mu.diff[1,1,3] <-mu[1,1,3]-mu[1,1,1]
  mu.diff[1,2,2] <-mu[1,2,2]-mu[1,1,1]
  mu.diff[1,2,3] <-mu[1,2,3]-mu[1,1,1]
  
  mu.diff[2,1,2] <-mu[2,1,2]-mu[2,1,1]
  mu.diff[2,1,3] <-mu[2,1,3]-mu[2,1,1]
  mu.diff[2,2,2] <-mu[2,2,2]-mu[2,1,1]
  mu.diff[2,2,3] <-mu[2,2,3]-mu[2,1,1]
  
}    


#Model

model.BHM <- BHM_2cont.bugs

m.sd.prior<-c(2,1)
Prior.mu.ex <- array(sort(rep(c(0,10), times=6)), dim=c(2,3,2)) 
p.exch <- cbind(c(0,0.5,0.5), c(0,0.5,0.5))
Prior.rho.mu <- c(-1,1)
Prior.mu.nex <- array(rep(c(0,0), times=6), dim=c(3,2,2)) 
Prior.v.nex <- array(sort(rep(c(5,2), times=6),decreasing =TRUE), dim=c(3,2,2))
Prior.v.nex[1,1:2,1:2]<-sqrt(10)

Prior.sigma.1 <- c(0,10)
Prior.sigma.2 <- c(0,10)
Prior.rho.y <- c(-1,1)


MCMC  <- c(25000, 30000, 2, 2)

Nobs <- sum(n)
trt1 <- data.sim[1:6]$trt1
trt2 <- data.sim[1:6]$trt2

y.new <- cbind(data.sim[1:6]$`1`, 
               data.sim[1:6]$`2`,
               data.sim[1:6]$`3`,
               data.sim[1:6]$`4`)

data.BHM <- list("Nobs", "trt1", "trt2", "y.new",
                 "m.sd.prior", "Prior.mu.ex", "p.exch",
                 "Prior.rho.mu", "Prior.mu.nex", "Prior.v.nex", 
                 "Prior.sigma.1", "Prior.sigma.2", "Prior.rho.y")

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

param.BHM <- c("mu", "mu.diff","rho.y","sigma.1","sigma.2","rho.mu")

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

BHM.fit$BUGSoutput$summary



##########################################################################################  
#########################################################################################  
##(BJSM with  bivariate model for parameter and likelihood)
##########################################################################################

BJSM_2cont.bugs<- function(){
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
  }
  
  ##############################################################################
  #Likelihood variance-covariance 
  ##############################################################################
  for(s in 1:2){
    for(k in 1:3){
      
      sigma.lkd[s,k,1,1] <- pow(sigma.1[k],2)
      sigma.lkd[s,k,2,2] <- pow(sigma.2[k],2)
      sigma.lkd[s,k,1,2] <- rho.y[k,1]*sqrt(sigma.lkd[s,k,1,1])*sqrt(sigma.lkd[s,k,2,2])
      sigma.lkd[s,k,2,1] <- sigma.lkd[s,k,1,2]
      prec.lkd[s,k,1:2,1:2] <- inverse(sigma.lkd[s,k,,])
    }
  }
  
  for(k in 1:3){
    
    rho.y[k,1]~ dunif(-1, 1)
  }
  
  ##############################################################################
  #Likelihood calculation
  ##############################################################################
  for (i in 1:Nobs) {
    
    mu.y1[i,1] <- mu.vec[1,1,trt1[i]]
    mu.y1[i,2] <- mu.vec[2,1,trt1[i]]
    mu.y2[i,1] <- mu.vec[1,2,trt2[i]]
    mu.y2[i,2] <- mu.vec[2,2,trt2[i]]
    
    y.new[i,1:2] ~ dmnorm(mu.y1[i,1:2], prec.lkd[1,trt1[i],,])
    y.new[i,3:4] ~ dmnorm(mu.y2[i,1:2], prec.lkd[2,trt2[i],,])
  }
  
  #Adding global inference
  mu.diff[1,1,2] <-mu.vec[1,1,2]-mu.vec[1,1,1]
  mu.diff[1,1,3] <-mu.vec[1,1,3]-mu.vec[1,1,1]
  mu.diff[1,2,2] <-mu.vec[1,2,2]-mu.vec[1,1,1]
  mu.diff[1,2,3] <-mu.vec[1,2,3]-mu.vec[1,1,1]
  
  mu.diff[2,1,2] <-mu.vec[2,1,2]-mu.vec[2,1,1]
  mu.diff[2,1,3] <-mu.vec[2,1,3]-mu.vec[2,1,1]
  mu.diff[2,2,2] <-mu.vec[2,2,2]-mu.vec[2,1,1]
  mu.diff[2,2,3] <-mu.vec[2,2,3]-mu.vec[2,1,1]
  
}  

#Model

model.BJSM <- BJSM_2cont.bugs
data.BJSM <-  list("Nobs", "trt1", "trt2", "y.new", 
                   "Prior.sigma.1", "Prior.sigma.2","Prior.mu.ex")

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

param.BJSM <- c("mu", "mu.diff","rho.y","rho.mu","sigma.1","sigma.2","alpha","mu.vec")

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

BJSM.fit$BUGSoutput$summary

# Summary

#Calcualte rMSE:
calculate_rmse <- function(posterior_samples, true_values) {
  # Ensure posterior_samples is a matrix and true_values is a vector
  posterior_samples <- as.matrix(posterior_samples)
  true_values <- as.vector(true_values)
  
  # Number of parameters
  n_params <- ncol(posterior_samples)
  
  # Initialize RMSE vector
  rmse_values <- numeric(n_params)
  
  # Calculate RMSE for each parameter
  for (i in 1:n_params) {
    # Calculate squared errors
    squared_errors <- (posterior_samples[, i] - true_values[i])^2
    # Calculate mean of squared errors
    mean_squared_error <- mean(squared_errors)
    # Calculate RMSE
    rmse_values[i] <- sqrt(mean_squared_error)
  }
  
  return(rmse_values)
}

#-------------------------------------------------------------------------------
#BHM
#-------------------------------------------------------------------------------

{
  # 
  # library(ggmcmc)
  # library(gridExtra)
  # 
  # fit.mcmc<-as.mcmc(BHM.fit)
  # 
  # S <- ggs(fit.mcmc)
  # trace_plots <- list(
  #   ggs_traceplot(S[S$Parameter == "mu[1,1,1]",]),
  #   ggs_traceplot(S[S$Parameter == "mu[1,1,2]",]),
  #   ggs_traceplot(S[S$Parameter == "mu[1,1,3]",]),
  #   ggs_traceplot(S[S$Parameter == "mu[1,2,2]",]),
  #   ggs_traceplot(S[S$Parameter == "mu[1,2,3]",]),
  #   ggs_traceplot(S[S$Parameter == "mu[2,1,1]",]),
  #   ggs_traceplot(S[S$Parameter == "mu[2,1,2]",]),
  #   ggs_traceplot(S[S$Parameter == "mu[2,1,3]",]),
  #   ggs_traceplot(S[S$Parameter == "mu[2,2,2]",]),
  #   ggs_traceplot(S[S$Parameter == "mu[2,2,3]",])
  # )
  # 
  # # Combine and save trace plots
  # combined_trace_plot <- marrangeGrob(trace_plots, nrow = 4, ncol = 3)
  # combined_trace_plot
  # 
  # modify_autocorr_plot <- function(plot, parameter) {
  #   plot +
  #     labs(subtitle = paste("Autocorrelation for", parameter)) +
  #     coord_cartesian(ylim = c(0, 1))  # Keep y-axis from 0 to 1
  # }
  # 
  # # List of parameters
  # parameters <- c(
  #   "mu[1,1,1]", "mu[1,1,2]", "mu[1,1,3]",
  #   "mu[1,2,2]", "mu[1,2,3]",
  #   "mu[2,1,1]", "mu[2,1,2]", "mu[2,1,3]",
  #   "mu[2,2,2]", "mu[2,2,3]"
  # )
  # 
  # # Create and modify autocorrelation plots
  # autocorr_plots <- lapply(parameters, function(param) {
  #   original_plot <- ggs_autocorrelation(S[S$Parameter == param, ])
  #   modify_autocorr_plot(original_plot, param)
  # })
  # 
  # # Combine and save autocorrelation plots
  # combined_autocorr_plot <- marrangeGrob(autocorr_plots, nrow = 4, ncol = 3)
  # combined_autocorr_plot
  
  #Output
  #--------------------------------------------------------------------------------------------------------------
  round(BHM.fit$BUGSoutput$summary,2)
  BHM=BHM.fit$BUGSoutput$sims.matrix
  BHM=cbind(BHM[,"mu.diff[1,1,2]"],BHM[,"mu.diff[1,1,3]"],BHM[,"mu.diff[2,1,2]"],BHM[,"mu.diff[2,1,3]"],
            BHM[,"mu[1,1,1]"], #outcome1 stage 1 p,l,h, stage 2 l h
            BHM[,"mu[2,1,1]"]) 
  BHM.mean=colMeans(BHM)
  #Bias
  BHM.Bias <- round(BHM.mean-mean.true,2)
  
  #Coverage:
  BHM.lower <- apply(BHM, 2, quantile, probs = 0.025)
  BHM.upper <- apply(BHM, 2, quantile, probs = 0.975)
  
  # Check if true value is inside the interval
  BHM.coverage <- (mean.true >= BHM.lower) & (mean.true <= BHM.upper)
  
  
  #Global inference:
  BHM.mu.diff=BHM.fit$BUGSoutput$sims.matrix[,c("mu.diff[1,1,2]",
                                                "mu.diff[1,1,3]",
                                                "mu.diff[2,1,2]",
                                                "mu.diff[2,1,3]")]
  
  ###Either low or high doses have significant different than placebo in either stages:
  ##Outcome 1
  BHM.out1.inf=(sum(BHM.mu.diff[,"mu.diff[1,1,2]"]>0 | BHM.mu.diff[,"mu.diff[1,1,3]"]>0)/nrow(BHM.mu.diff))>0.975
  BHM.out1.l.inf=(sum(BHM.mu.diff[,"mu.diff[1,1,2]"]>0)/nrow(BHM.mu.diff))>0.975
  BHM.out1.h.inf=(sum(BHM.mu.diff[,"mu.diff[1,1,3]"]>0)/nrow(BHM.mu.diff))>0.975
  ##Outcome 2
  BHM.out2.inf=(sum(BHM.mu.diff[,"mu.diff[2,1,2]"]>0 | BHM.mu.diff[,"mu.diff[2,1,3]"]>0)/nrow(BHM.mu.diff))>0.975
  BHM.out2.l.inf=(sum(BHM.mu.diff[,"mu.diff[2,1,2]"]>0)/nrow(BHM.mu.diff))>0.975
  BHM.out2.h.inf=(sum(BHM.mu.diff[,"mu.diff[2,1,3]"]>0)/nrow(BHM.mu.diff))>0.975
  ###Either low or high doses have significant different than placebo in either outcomes stages:
  BHM.all.inf=(sum((BHM.mu.diff[,"mu.diff[1,1,2]"]>0 & 
                      BHM.mu.diff[,"mu.diff[2,1,2]"]>0) | 
                     (BHM.mu.diff[,"mu.diff[1,1,3]"]>0 & 
                        BHM.mu.diff[,"mu.diff[2,1,3]"]>0)  )/nrow(BHM.mu.diff))>0.975
  
  BHM.all.l.inf=(sum(BHM.mu.diff[,"mu.diff[1,1,2]"]>0 & 
                       BHM.mu.diff[,"mu.diff[2,1,2]"]>0 )/nrow(BHM.mu.diff))>0.975
  
  
  BHM.all.h.inf=(sum(BHM.mu.diff[,"mu.diff[1,1,3]"]>0 & 
                       BHM.mu.diff[,"mu.diff[2,1,3]"]>0  )/nrow(BHM.mu.diff))>0.975
  BHM.either.inf=(sum(BHM.mu.diff[,"mu.diff[1,1,2]"]>0 | 
                        BHM.mu.diff[,"mu.diff[1,1,3]"]>0 | 
                        BHM.mu.diff[,"mu.diff[2,1,2]"]>0 | 
                        BHM.mu.diff[,"mu.diff[2,1,3]"]>0  )/nrow(BHM.mu.diff))>0.975
  BHM.either.l.inf=(sum(BHM.mu.diff[,"mu.diff[1,1,2]"]>0 | 
                          BHM.mu.diff[,"mu.diff[2,1,2]"]>0 )/nrow(BHM.mu.diff))>0.975
  BHM.either.h.inf=(sum(BHM.mu.diff[,"mu.diff[1,1,3]"]>0 | 
                          BHM.mu.diff[,"mu.diff[2,1,3]"]>0  )/nrow(BHM.mu.diff))>0.975
  BHM.inf=c(BHM.out1.l.inf,BHM.out1.h.inf,BHM.out1.inf,
            BHM.out2.l.inf,BHM.out2.h.inf,BHM.out2.inf,
            BHM.all.l.inf,BHM.all.h.inf,BHM.all.inf,
            BHM.either.l.inf,BHM.either.h.inf,BHM.either.inf)
  BHM.absbias=abs(BHM.Bias)
}

#-------------------------------------------------------------------------------
#BJSM
#-------------------------------------------------------------------------------
{
  BJSM.fit$BUGSoutput$summary
  BJSM=BJSM.fit$BUGSoutput$sims.matrix
  BJSM=cbind(BJSM[,"mu.diff[1,1,2]"],BJSM[,"mu.diff[1,1,3]"],BJSM[,"mu.diff[2,1,2]"],BJSM[,"mu.diff[2,1,3]"],
             BJSM[,"mu[1,1]"], #outcome1 stage 1 p,l,h, stage 2 l h
             BJSM[,"mu[2,1]"]) #outcome2 stage 1 p,l,h, stage 2 l h
  BJSM.mean=colMeans(BJSM)
  
  #True mean (only use sta)
  BJSM.Bias <- round(BJSM.mean-mean.true,2)
  #Coverage:
  BJSM.lower <- apply(BJSM, 2, quantile, probs = 0.025)
  BJSM.upper <- apply(BJSM, 2, quantile, probs = 0.975)
  
  # Check if true value is inside the interval
  BJSM.coverage <- (mean.true >= BJSM.lower) & (mean.true <= BJSM.upper)
  #Global inference:
  BJSM.mu.diff=BJSM.fit$BUGSoutput$sims.matrix[,c("mu.diff[1,1,2]",
                                                  "mu.diff[1,1,3]",
                                                  "mu.diff[2,1,2]",
                                                  "mu.diff[2,1,3]")]
  
  ###Either low or high doses have significant different than placebo in either stages:
  ##Outcome 1
  BJSM.out1.inf=(sum(BJSM.mu.diff[,"mu.diff[1,1,2]"]>0 | BJSM.mu.diff[,"mu.diff[1,1,3]"]>0)/nrow(BJSM.mu.diff))>0.975
  BJSM.out1.l.inf=(sum(BJSM.mu.diff[,"mu.diff[1,1,2]"]>0)/nrow(BJSM.mu.diff))>0.975
  BJSM.out1.h.inf=(sum(BJSM.mu.diff[,"mu.diff[1,1,3]"]>0)/nrow(BJSM.mu.diff))>0.975
  ##Outcome 2
  BJSM.out2.inf=(sum(BJSM.mu.diff[,"mu.diff[2,1,2]"]>0 | BJSM.mu.diff[,"mu.diff[2,1,3]"]>0)/nrow(BJSM.mu.diff))>0.975
  BJSM.out2.l.inf=(sum(BJSM.mu.diff[,"mu.diff[2,1,2]"]>0)/nrow(BJSM.mu.diff))>0.975
  BJSM.out2.h.inf=(sum(BJSM.mu.diff[,"mu.diff[2,1,3]"]>0)/nrow(BJSM.mu.diff))>0.975
  ###Either low or high doses have significant different than placebo in either outcomes stages:
  BJSM.all.inf=(sum((BJSM.mu.diff[,"mu.diff[1,1,2]"]>0 & 
                       BJSM.mu.diff[,"mu.diff[2,1,2]"]>0) | 
                      (BJSM.mu.diff[,"mu.diff[1,1,3]"]>0 & 
                         BJSM.mu.diff[,"mu.diff[2,1,3]"]>0)  )/nrow(BJSM.mu.diff))>0.975
  
  BJSM.all.l.inf=(sum(BJSM.mu.diff[,"mu.diff[1,1,2]"]>0 & 
                        BJSM.mu.diff[,"mu.diff[2,1,2]"]>0 )/nrow(BJSM.mu.diff))>0.975
  
  
  BJSM.all.h.inf=(sum(BJSM.mu.diff[,"mu.diff[1,1,3]"]>0 & 
                        BJSM.mu.diff[,"mu.diff[2,1,3]"]>0  )/nrow(BJSM.mu.diff))>0.975
  BJSM.either.inf=(sum(BJSM.mu.diff[,"mu.diff[1,1,2]"]>0 | 
                         BJSM.mu.diff[,"mu.diff[1,1,3]"]>0 | 
                         BJSM.mu.diff[,"mu.diff[2,1,2]"]>0 | 
                         BJSM.mu.diff[,"mu.diff[2,1,3]"]>0  )/nrow(BJSM.mu.diff))>0.975
  BJSM.either.l.inf=(sum(BJSM.mu.diff[,"mu.diff[1,1,2]"]>0 | 
                           BJSM.mu.diff[,"mu.diff[2,1,2]"]>0 )/nrow(BJSM.mu.diff))>0.975
  BJSM.either.h.inf=(sum(BJSM.mu.diff[,"mu.diff[1,1,3]"]>0 | 
                           BJSM.mu.diff[,"mu.diff[2,1,3]"]>0  )/nrow(BJSM.mu.diff))>0.975
  BJSM.inf=c(BJSM.out1.l.inf,BJSM.out1.h.inf,BJSM.out1.inf,
             BJSM.out2.l.inf,BJSM.out2.h.inf,BJSM.out2.inf,
             BJSM.all.l.inf,BJSM.all.h.inf,BJSM.all.inf,
             BJSM.either.l.inf,BJSM.either.h.inf,BJSM.either.inf)
  BJSM.absbias=abs(BJSM.Bias)
}