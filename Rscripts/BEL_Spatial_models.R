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
require(rgdal)          # For readOGR()
require(rgeos)          # For unionSpatialPolygons(), gIntersection(), gBuffer()
require(maptools)       # For unionSpatialPolygons()
require(ggplot2)        # For fortify(), ggplot()
require(readxl)         # For read_excel()
require(magrittr)       # For the pipe operator %>%
require(scales)         # For rescale()
require(dplyr)          # For inner_join(), bind_rows(),### between(), mutate()
require(gridExtra)      # For grid.arrange()
require(tidyr)          # For gather()
require(spdep) # for reading .gal files)
require(MRH) #

map <- readOGR("Data/scotlip/scotlip.shp", verbose = FALSE)

# Create border map from full map
N <- length(map)
map.border <- unionSpatialPolygons(map, IDs = rep(1, N))
map.border <- fortify(map.border)

# Fortify
map.df <- fortify(map)

# Make shapefile dataframe IDs numeric
map.df$id <- as.numeric(map.df$id)
if(!all(range(map.df$id) == c(1, N))){
  map.df$id <- as.numeric(map.df$id) + 1  # Need to add 1 so range is 1:N
  print(range(map.df$id))
}
stopifnot(all(unique(map.df$id) == 1:N))

# reading data
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
#save(W,file="Data/scotlip/W.RDS")
load("Data/scotlip/W.RDS")

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

# fitting BEL BYM model taking rho= 1
library(parallel)
cluster<-makeCluster(3)
clusterEvalQ(cl=cluster,.libPaths("c:/software/Rpackages"))
clusterEvalQ(cl=cluster,library(BELSpatial))
clusterExport(cl=cluster,varlist = c("y","x","n","p","var","beta_init", "psi_init", "tau_init"
                                     ,"R", "wi"))
BEL_BYM_scotlip<-clusterApply(cl=cluster, x=1:3, function(z){BEL_leroux_new(y,x,n,p,var,rho=1,niter=500,
                                                          beta_init, psi_init, tau_init,R, wi, sd_psi=0.35, 
                                                           sd_beta=1, sd_tau=0.4)})
save(BEL_BYM_scotlip,file="Results/BEL_BYM_scotlip.RData")

# fitting BEL BYM model taking rho= 0.75

BEL_leroux_scotlip<-clusterApply(cl=cluster, x=1:3, function(z){BEL_leroux_new(y,x,n,p,var,rho=0.75,niter=500,
                                                                            beta_init, psi_init, tau_init,R, wi, sd_psi=0.35, 
                                                                            sd_beta=1, sd_tau=0.4)})
save(BEL_leroux_scotlip,file="Results/BEL_Leroux_scotlip.RData")

# Porter's BSHEL model
B<-W
B_plus<-diag(rowSums(B))
M=M_create(y,x,B)
MBM=MBM_create(M,B,B_plus)
q=dim(MBM)[2]
psi_init <- rep(0,q) 
wi=wi_init
var<- as.numeric(var(y- x%*%beta_init))
beta_mele<- mele( x = x, tet= beta_init,y=y,var=var) 
mu_init<- x%*% beta_mele + M%*%psi_init
beta_init<-beta_mele

wi_mu<- el.test(y-mu_init,0)$wts # computing el weights using emplik package
wi_mu<-wi_mu/sum(wi_mu) # sum(wi) = 1 and wi>0 constraints 
wi<-wi_mu

#Fitting the Porter BSHEL model

clusterExport(cl=cluster,varlist = c("y","x","n","p","var","beta_init", "psi_init", "tau_init"
                                     ,"B","B_plus","q","M","MBM", "wi"))
Porter_BSHEL_scotlip<-clusterApply(cl=cluster, x=1:3, fun= function(z){BSHEL(y,x,n,p,q,var,niter=500,beta_init, 
                                                                     psi_init, tau_init,M,MBM, wi, 
                                                                     sd_psi=0.00001, 
                                                                     sd_beta=0.0001, sd_tau=0.9)})

save(Porter_BSHEL_scotlip,file="Results/BSHEL_porter_scotlip.RData")
stopCluster(cl=cluster)
