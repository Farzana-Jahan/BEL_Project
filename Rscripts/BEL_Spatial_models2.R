# for office PC only run the /libpaths code, otherwise comment it out
# .libPaths("c:/software/Rpackages")
# devtools::install_github("danwkenn/BELSpatial") # installing the package
#devtools::install_local("C:\\R dir\\Leroux-Empirical-Likelihood-master\\Leroux-Empirical-Likelihood-master\\BELSpatial",
#force=TRUE)
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
## Loading the North Carolina SIDS Data
map2 <- readOGR(system.file("shapes/sids.shp", package = "spData")[1])
fp.wd<- getwd()
# Add a projection
proj4string(map2) <- CRS("+proj=longlat +datum=NAD27")

# Fortify for plotting
map2.df <- fortify(map2)

# Make shapefile dataframe IDs numeric
N <- length(map2)
map.df$id <- as.numeric(map2.df$id)
if(!all(range(map2.df$id) == c(1, N))){
  map2.df$id <- as.numeric(map2.df$id) + 1  # Need to add 1 so range is 1:N
}
stopifnot(all(unique(map2.df$id) == 1:N))
# Create border map from full map
map2.border <- unionSpatialPolygons(map2, IDs = rep(1, N))
# View data contained within shapefile
data <- as.data.frame(map2)
data <- data[,-c(1:7, 15:22)]
dat.1 <- data.frame(data[,c(1:4)], year = "1974-78")
names(dat.1) <- c("area.ID", "pop", "observed", "x", "year")
dat.2 <- data.frame(data[,c(1, 5:7)], year = "1979-84")
names(dat.2) <- names(dat.1)
data <- rbind(dat.1, dat.2)
rm(dat.1, dat.2)
data$expected <- c(
  data$pop[1:100] * sum(data$observed[1:100]) / sum(data$pop[1:100]),
  data$pop[101:200] * sum(data$observed[101:200]) / sum(data$pop[101:200])
)
data$raw <- data$observed / data$expected
data_1<-data[data$year=="1974-78",]
# creating neighbourhood matrix
sids_nb<-read.gal("Data/SIDS Data/shapefile_SIDS.gal", override.id = TRUE)
#class(sids_nb)
W<-nb2mat(sids_nb,style="B")
nblist<-nb2listw(sids_nb)
#creating symmetric neighbourhood matrix for BYM in CARBAYES
rownames(W)<-c()
ind <- upper.tri(W)
W[ind] <- t(W)[ind] 
ni<-rowSums(W) # no. of neighbours for each area
R2<-diag(ni)
for(i in 1:nrow(R2))
{
  R2[i,which(W[i,]==1)]<- -1
}
R<- R2
# fitting BEl spatial model utilising the scaled/standardised covariate
x<- cbind(1, scale(data_1$x))

# using log(rawSIRs) as the response
data_1$raw[data_1$raw==0]<-0.1
y= log(data_1$raw)
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
beta<-beta_mele
wi_mu<- el.test(y-mu_init,0)$wts # computing el weights using emplik package
wi_mu<-wi_mu/sum(wi_mu) # sum(wi) = 1 and wi>0 constraints 
wi<-wi_mu

# fitting BEL BYM model takign rho= 1
BEL_BYM_chain1_SIDS<- BEL_leroux_new(y,x,n,p,var,rho=1,niter=10000,beta_init, psi_init, tau_init,R, wi, 
                                     sd_psi=0.002, sd_beta=0.009, sd_tau=4) # how to parallelise this code? or run multple codes with different initial values and save to get the posteriors?
save(BEL_BYM_chain1_SIDS,file="Results/BYM_BEL1_SIDS.RData")
# fitting BEL leroux model takign rho= 1
BEL_leroux_chain1_SIDS<- BEL_leroux_new(y,x,n,p,var,rho=0.75,niter=5000000,beta_init, psi_init, tau_init,R, wi, sd_psi=100, 
                                        sd_beta=100, sd_tau=4) # how to parallelise this code? or run multple codes with different initial values and save to get the posteriors?
save(BEL_leroux_chain1_SIDS,file="Results/Leroux_BEL1_SIDS.RData")
# Porter's BSHEL model
B<-W
B_plus<-diag(rowSums(B))
M=M_create(y,x,B)
MBM=MBM_create(M,B,B_plus)
q=dim(MBM)[2]
psi_init <- rep(0,q) 
wi=wi_init

beta_mele<- mele( x = x, tet= beta_init,y=y,var=var) 
mu_init<- x%*% beta_mele + M%*%psi_init
beta_init<-beta_mele
var<- as.numeric(var(y- x%*%beta_init))
wi_mu<- el.test(y-mu_init,0)$wts # computing el weights using emplik package
wi_mu<-wi_mu/sum(wi_mu) # sum(wi) = 1 and wi>0 constraints 
wi<-wi_mu

#Fitting the POrter BSHEL model
porter_chain1_SIDS<-BSHEL(y,x,n,p,var,niter=10000,beta_init, 
                          psi_init, tau_init,M,MBM, wi, sd_psi=0.00001, sd_beta=0.00001, sd_tau=0.9)
save(porter_chain1_SIDS,file="Results/BSHEL_porter_SIDS.RData")

