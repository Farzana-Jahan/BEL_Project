install.packages("parallel")
library(parallel)

cluster<-makeCluster(3)
clusterEvalQ(cl=cluster,.libPaths("c:/software/Rpackages"))
clusterEvalQ(cl=cluster,library(BELSpatial))
clusterExport(cl=cluster,varlist = c("y","x","n","p","var",rho=1,niter=5,"beta_init", "psi_init", "tau_init"
                                     ,"R", "wi"))
clusterApply(cl=cluster, x= list(1,2,3), fun=function(x){x^2})
stopCluster(cl=cluster)

# qsub jobs/generic.job.sh
qstat -u n9863176