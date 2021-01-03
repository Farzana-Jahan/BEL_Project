# for office PC only run the /libpaths code, otherwise comment it out
# .libPaths("c:/software/Rpackages")
#devtools::install_github("danwkenn/BELSpatial") # installing the package
#devtools::install_local("C:\\R dir\\Leroux-Empirical-Likelihood-master\\Leroux-Empirical-Likelihood-master\\BELSpatial",force=TRUE)
# loading the data (Scottish Lip Cancer Data)

# libraries
require(BELSpatial)
require(gmm)
require(emplik)
require(MASS)
require(ggplot2)        # For fortify(), ggplot()
require(readxl)         # For read_excel()
require(magrittr)       # For the pipe operator %>%
require(scales)         # For rescale()
require(dplyr)          # For inner_join(), bind_rows(),### between(), mutate()
require(gridExtra)      # For grid.arrange()
require(tidyr)          # For gather()
library(tidyverse)



# reading excel data
data1 <- read_excel("Data/scotlip/scotlip.xlsx") %>% 
  data.frame %>%
  subset(select = -c(1:4, 6:7, 9))

# reading neighbourhood matrix from text file

#scot_nb<-read.gal("Data/scotlip/scotlip.gal", override.id = TRUE)
#class(scot_nb)
#W<-nb2mat(scot_nb,style="B")
#nblist<-nb2listw(scot_nb)
#creating symmetric neighbourhood matrix for BYM in CARBAYES
#rownames(W)<-c()
#ind <- upper.tri(W)
#W[ind] <- t(W)[ind] 
#saveRDS(W,"Data/scotlip/W.RDS")
W<-readRDS("Data/scotlip/W.RDS")

ni<-rowSums(W) # no. of neighbours for each area
R<-diag(ni)
for(i in 1:nrow(R))
{
  R[i,which(W[i,]==1)]<- -1
}

# calculating SIR and adding a column to data
data1$CANCER[data1$CANCER==0]<- 0.1
data1$SIR<- data1$CANCER/data1$CEXP

# fitting BEl spatial model utilising the scaled/standardised covariate
x<- cbind(1, scale(data1$AFF))

# using log(rawSIRs) as the response
y= log(data1$SIR)
# initial values needed before fitting models

n<- length(y) # no. of observations
p<- dim(x)[2] # no. of covariates
alpha_1<-1 # hyperparamter for tau prior
alpha_2<-0.01 # hyperparamter for tau prior
tau_inv_init<- rgamma(1,alpha_1,alpha_2) # using IG prior(1,1) for tau_inv
tau_init<- 1/tau_inv_init
g<- 10# G prior evaluated at 10 for regression coefficients' prior (Zellner prior)
prior_mean_beta<- rep(0,p) # p is the number of regression parameters, in case of one covariate, p=2
beta_init<- rnorm(2,prior_mean_beta, (1/g)*tau_inv_init)
wi_init<- 1/length(y) # y be the response variable from the data
psi_init <- rep(0,n)
var<- as.numeric(var(y- x%*%beta_init))
# calculating MELE of Beta, beta_mele
wi=wi_init
beta_mele<- mele(x,tet=beta_init,y=y,var=var)
mu_init<- x%*% beta_mele + psi_init
beta_init<-beta_mele
wi_mu<- el.test(y-mu_init,0)$wts # computing el weights using emplik package
wi_mu<-wi_mu/sum(wi_mu) # sum(wi) = 1 and wi>0 constraints 
wi<-wi_mu


# fitting BEL BYM model taking rho= 0.75

BEL_leroux_scotlip<-BEL_leroux_new(y,x,n,p,var,rho=0.1,niter=10000,
                                  beta_init, psi_init, tau_init,R, wi, sd_psi=0.35, 
                                  sd_beta=1, sd_tau=0.4)
Beta_leroux_BEL<- rbind(BEL_leroux_scotlip$Beta[1,],BEL_leroux_scotlip$Beta[1,])
psi_leroux_BEL<-matrix(c(BEL_leroux_scotlip$psi[,1]),nrow=10000, ncol=1)
for(i in 2:56){
  psi_leroux_BEL<-cbind(psi_leroux_BEL,BEL_leroux_scotlip$psi[,i])
}
tau_leroux_BEL<-BEL_leroux_scotlip$tau
theta<-t(Beta_leroux_BEL)%*%t(x)+psi_leroux_BEL
source("Rscripts/waic.R")
WAIC_leroux_BEL_0.1<-get.WAIC.BEL(theta=theta,y=y,x=x)
BEL_leroux_scotlip<-BEL_leroux_new(y,x,n,p,var,rho=0.2,niter=10000,
                                   beta_init, psi_init, tau_init,R, wi, sd_psi=0.35, 
                                   sd_beta=1, sd_tau=0.4)
Beta_leroux_BEL<- rbind(BEL_leroux_scotlip$Beta[1,],BEL_leroux_scotlip$Beta[1,])
psi_leroux_BEL<-matrix(c(BEL_leroux_scotlip$psi[,1]),nrow=10000, ncol=1)
for(i in 2:56){
  psi_leroux_BEL<-cbind(psi_leroux_BEL,BEL_leroux_scotlip$psi[,i])
}
tau_leroux_BEL<-BEL_leroux_scotlip$tau
theta<-t(Beta_leroux_BEL)%*%t(x)+psi_leroux_BEL
source("Rscripts/waic.R")
WAIC_leroux_BEL_0.2<-get.WAIC.BEL(theta=theta,y=y,x=x)
BEL_leroux_scotlip<-BEL_leroux_new(y,x,n,p,var,rho=0.3,niter=10000,
                                   beta_init, psi_init, tau_init,R, wi, sd_psi=0.35, 
                                   sd_beta=1, sd_tau=0.4)
Beta_leroux_BEL<- rbind(BEL_leroux_scotlip$Beta[1,],BEL_leroux_scotlip$Beta[1,])
psi_leroux_BEL<-matrix(c(BEL_leroux_scotlip$psi[,1]),nrow=10000, ncol=1)
for(i in 2:56){
  psi_leroux_BEL<-cbind(psi_leroux_BEL,BEL_leroux_scotlip$psi[,i])
}
tau_leroux_BEL<-BEL_leroux_scotlip$tau
theta<-t(Beta_leroux_BEL)%*%t(x)+psi_leroux_BEL
source("Rscripts/waic.R")
WAIC_leroux_BEL_0.3<-get.WAIC.BEL(theta=theta,y=y,x=x)
BEL_leroux_scotlip<-BEL_leroux_new(y,x,n,p,var,rho=0.4,niter=10000,
                                   beta_init, psi_init, tau_init,R, wi, sd_psi=0.35, 
                                   sd_beta=1, sd_tau=0.4)
Beta_leroux_BEL<- rbind(BEL_leroux_scotlip$Beta[1,],BEL_leroux_scotlip$Beta[1,])
psi_leroux_BEL<-matrix(c(BEL_leroux_scotlip$psi[,1]),nrow=10000, ncol=1)
for(i in 2:56){
  psi_leroux_BEL<-cbind(psi_leroux_BEL,BEL_leroux_scotlip$psi[,i])
}
tau_leroux_BEL<-BEL_leroux_scotlip$tau
theta<-t(Beta_leroux_BEL)%*%t(x)+psi_leroux_BEL
source("Rscripts/waic.R")
WAIC_leroux_BEL_0.4<-get.WAIC.BEL(theta=theta,y=y,x=x)
BEL_leroux_scotlip<-BEL_leroux_new(y,x,n,p,var,rho=0.5,niter=10000,
                                   beta_init, psi_init, tau_init,R, wi, sd_psi=0.35, 
                                   sd_beta=1, sd_tau=0.4)
Beta_leroux_BEL<- rbind(BEL_leroux_scotlip$Beta[1,],BEL_leroux_scotlip$Beta[1,])
psi_leroux_BEL<-matrix(c(BEL_leroux_scotlip$psi[,1]),nrow=10000, ncol=1)
for(i in 2:56){
  psi_leroux_BEL<-cbind(psi_leroux_BEL,BEL_leroux_scotlip$psi[,i])
}
tau_leroux_BEL<-BEL_leroux_scotlip$tau
theta<-t(Beta_leroux_BEL)%*%t(x)+psi_leroux_BEL
source("Rscripts/waic.R")
WAIC_leroux_BEL_0.5<-get.WAIC.BEL(theta=theta,y=y,x=x)
BEL_leroux_scotlip<-BEL_leroux_new(y,x,n,p,var,rho=0.6,niter=10000,
                                   beta_init, psi_init, tau_init,R, wi, sd_psi=0.35, 
                                   sd_beta=1, sd_tau=0.4)
Beta_leroux_BEL<- rbind(BEL_leroux_scotlip$Beta[1,],BEL_leroux_scotlip$Beta[1,])
psi_leroux_BEL<-matrix(c(BEL_leroux_scotlip$psi[,1]),nrow=10000, ncol=1)
for(i in 2:56){
  psi_leroux_BEL<-cbind(psi_leroux_BEL,BEL_leroux_scotlip$psi[,i])
}
tau_leroux_BEL<-BEL_leroux_scotlip$tau
theta<-t(Beta_leroux_BEL)%*%t(x)+psi_leroux_BEL
source("Rscripts/waic.R")
WAIC_leroux_BEL_0.6<-get.WAIC.BEL(theta=theta,y=y,x=x)
BEL_leroux_scotlip<-BEL_leroux_new(y,x,n,p,var,rho=0.7,niter=10000,
                                   beta_init, psi_init, tau_init,R, wi, sd_psi=0.35, 
                                   sd_beta=1, sd_tau=0.4)
Beta_leroux_BEL<- rbind(BEL_leroux_scotlip$Beta[1,],BEL_leroux_scotlip$Beta[1,])
psi_leroux_BEL<-matrix(c(BEL_leroux_scotlip$psi[,1]),nrow=10000, ncol=1)
for(i in 2:56){
  psi_leroux_BEL<-cbind(psi_leroux_BEL,BEL_leroux_scotlip$psi[,i])
}
tau_leroux_BEL<-BEL_leroux_scotlip$tau
theta<-t(Beta_leroux_BEL)%*%t(x)+psi_leroux_BEL
source("Rscripts/waic.R")
WAIC_leroux_BEL_0.7<-get.WAIC.BEL(theta=theta,y=y,x=x)
BEL_leroux_scotlip<-BEL_leroux_new(y,x,n,p,var,rho=0.8,niter=10000,
                                   beta_init, psi_init, tau_init,R, wi, sd_psi=0.35, 
                                   sd_beta=1, sd_tau=0.4)
Beta_leroux_BEL<- rbind(BEL_leroux_scotlip$Beta[1,],BEL_leroux_scotlip$Beta[1,])
psi_leroux_BEL<-matrix(c(BEL_leroux_scotlip$psi[,1]),nrow=10000, ncol=1)
for(i in 2:56){
  psi_leroux_BEL<-cbind(psi_leroux_BEL,BEL_leroux_scotlip$psi[,i])
}
tau_leroux_BEL<-BEL_leroux_scotlip$tau
theta<-t(Beta_leroux_BEL)%*%t(x)+psi_leroux_BEL
source("Rscripts/waic.R")
WAIC_leroux_BEL_0.8<-get.WAIC.BEL(theta=theta,y=y,x=x)
BEL_leroux_scotlip<-BEL_leroux_new(y,x,n,p,var,rho=0.9,niter=10000,
                                   beta_init, psi_init, tau_init,R, wi, sd_psi=0.35, 
                                   sd_beta=1, sd_tau=0.4)
Beta_leroux_BEL<- rbind(BEL_leroux_scotlip$Beta[1,],BEL_leroux_scotlip$Beta[1,])
psi_leroux_BEL<-matrix(c(BEL_leroux_scotlip$psi[,1]),nrow=10000, ncol=1)
for(i in 2:56){
  psi_leroux_BEL<-cbind(psi_leroux_BEL,BEL_leroux_scotlip$psi[,i])
}
tau_leroux_BEL<-BEL_leroux_scotlip$tau
theta<-t(Beta_leroux_BEL)%*%t(x)+psi_leroux_BEL
source("Rscripts/waic.R")
WAIC_leroux_BEL_0.9<-get.WAIC.BEL(theta=theta,y=y,x=x)
WAIC_scotlip<- data.frame(rho=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9), WAIC=c(WAIC_leroux_BEL_0.1,WAIC_leroux_BEL_0.2,
                                                                             WAIC_leroux_BEL_0.3,WAIC_leroux_BEL_0.4,
                                                                             WAIC_leroux_BEL_0.5,WAIC_leroux_BEL_0.6,
                                                                             WAIC_leroux_BEL_0.7,WAIC_leroux_BEL_0.8,
                                                                             WAIC_leroux_BEL_0.9))
write_csv(WAIC_scotlip,"Results/WAIC_scotlip.csv")