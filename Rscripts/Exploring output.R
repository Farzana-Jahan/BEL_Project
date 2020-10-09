# MCMC samples are recorded in Results folder
# to get trace plots for all three parallel chains
# BYM_BEL
load("Results/BEL_BYM_scotlip.RData")
plot(1:900000,BEL_BYM_scotlip[[1]]$Beta[1,100001:1000000],type="l",main="Intercept",col="blue")
lines(1:900000,BEL_BYM_scotlip[[2]]$Beta[1,100001:1000000],type="l",main="Intercept",col="red")
lines(1:900000,BEL_BYM_scotlip[[3]]$Beta[1,100001:1000000],type="l",main="Intercept",col="green")
plot(1:900000,BEL_BYM_scotlip[[1]]$Beta[2,100001:1000000],type="l",main="Intercept",col="blue")
lines(1:900000,BEL_BYM_scotlip[[2]]$Beta[2,100001:1000000],type="l",main="Intercept",col="red")
lines(1:900000,BEL_BYM_scotlip[[3]]$Beta[2,100001:1000000],type="l",main="Intercept",col="green")
plot(1:900000,BEL_BYM_scotlip[[1]]$tau[100001:1000000],type="l",main="Precision paramter",col="blue")
lines(1:900000,BEL_BYM_scotlip[[2]]$tau[100001:1000000],type="l",main="Precision paramter",col="red")
lines(1:900000,BEL_BYM_scotlip[[3]]$tau[100001:1000000],type="l",main="Precision paramter",col="green")
plot(1:900000,BEL_BYM_scotlip[[1]]$psi[1,100001:1000000],type="l",main="Spatial RE",col="blue")
lines(1:900000,BEL_BYM_scotlip[[2]]$psi[1,100001:1000000],type="l",main="Spatial RE",col="red")
lines(1:900000,BEL_BYM_scotlip[[3]]$psi[1,100001:1000000],type="l",main="Spatial RE",col="green")
# Leroux_BEL model
load("Results/BEL_Leroux_scotlip.RData")
plot(1:900000,BEL_Leroux_scotlip[[1]]$Beta[1,100001:1000000],type="l",main="Intercept",col="blue")
lines(1:900000,BEL_Leroux_scotlip[[2]]$Beta[1,100001:1000000],type="l",main="Intercept",col="red")
lines(1:900000,BEL_Leroux_scotlip[[3]]$Beta[1,100001:1000000],type="l",main="Intercept",col="green")
plot(1:900000,BEL_Leroux_scotlip[[1]]$Beta[2,100001:1000000],type="l",main="Intercept",col="blue")
lines(1:900000,BEL_Leroux_scotlip[[2]]$Beta[2,100001:1000000],type="l",main="Intercept",col="red")
lines(1:900000,BEL_Leroux_scotlip[[3]]$Beta[2,100001:1000000],type="l",main="Intercept",col="green")
plot(1:900000,BEL_Leroux_scotlip[[1]]$tau[100001:1000000],type="l",main="Precision paramter",col="blue")
lines(1:900000,BEL_Leroux_scotlip[[2]]$tau[100001:1000000],type="l",main="Precision paramter",col="red")
lines(1:900000,BEL_Leroux_scotlip[[3]]$tau[100001:1000000],type="l",main="Precision paramter",col="green")
plot(1:900000,BEL_Leroux_scotlip[[1]]$psi[1,100001:1000000],type="l",main="Spatial RE",col="blue")
lines(1:900000,BEL_Leroux_scotlip[[2]]$psi[1,100001:1000000],type="l",main="Spatial RE",col="red")
lines(1:900000,BEL_Leroux_scotlip[[3]]$psi[1,100001:1000000],type="l",main="Spatial RE",col="green")

# Porter_BSHEL
load("Results/BSHEL_porter_scotlip.RData")
plot(1:900000,Porter_BSHEL_scotlip[[1]]$Beta[1,100001:1000000],type="l",main="Intercept",col="blue")
lines(1:900000,Porter_BSHEL_scotlip[[2]]$Beta[1,100001:1000000],type="l",main="Intercept",col="red")
lines(1:900000,Porter_BSHEL_scotlip[[3]]$Beta[1,100001:1000000],type="l",main="Intercept",col="green")
plot(1:900000,Porter_BSHEL_scotlip[[1]]$Beta[2,100001:1000000],type="l",main="Intercept",col="blue")
lines(1:900000,Porter_BSHEL_scotlip[[2]]$Beta[2,100001:1000000],type="l",main="Intercept",col="red")
lines(1:900000,Porter_BSHEL_scotlip[[3]]$Beta[2,100001:1000000],type="l",main="Intercept",col="green")
plot(1:900000,Porter_BSHEL_scotlip[[1]]$tau[100001:1000000],type="l",main="Precision paramter",col="blue")
lines(1:900000,Porter_BSHEL_scotlip[[2]]$tau[100001:1000000],type="l",main="Precision paramter",col="red")
lines(1:900000,Porter_BSHEL_scotlip[[1]]$tau[100001:1000000],type="l",main="Precision paramter",col="green")
plot(1:900000,Porter_BSHEL_scotlip[[1]]$psi[1,100001:1000000],type="l",main="Spatial RE",col="blue")
lines(1:900000,Porter_BSHEL_scotlip[[2]]$psi[1,100001:1000000],type="l",main="Spatial RE",col="red")
lines(1:900000,Porter_BSHEL_scotlip[[1]]$psi[1,100001:1000000],type="l",main="Spatial RE",col="green")




