# User-defined function for computing WAIC#calculate theta as x%*%beta
library(emplik)
get.WAIC.BEL <- function(theta, y,x){
  # Pointwise evaluation of the likelihood and log-likelihood
  Like <- matrix(NA,nrow(theta), ncol(theta)) # ncol= n, observations, nrow= iterations
  for(i in 1:nrow(theta)){
    Like[i,]<-el.test(y-theta[i,],0)$wts/sum(el.test(y-theta[i,],0)$wts)
  }
  logLike <- log(Like)
  Var <- apply(logLike, 2, var)
  S <- apply(Like, 2, mean)
  WAIC <- 2 * sum(Var - log(S))
  return(WAIC)
}



# User-defined function for computing WAIC
# row is number of  iterations and col is number of observations here for theta
get.WAIC <- function(theta, y, L){
  # Pointwise evaluation of the likelihood and log-likelihood
  LogLike <- matrix(NA, nrow=nrow(theta), ncol=ncol(theta))
  for(i in 1:ncol(theta)){
    LogLike[,i] <- L(y[i], theta[,i], log = TRUE)
  }
  Like <- exp(LogLike)
  Var <- apply(LogLike, 2, var)
  S <- apply(Like, 2, mean)
  WAIC <- 2 * sum(Var - log(S))
  return(WAIC)
}
