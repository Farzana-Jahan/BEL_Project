# devtools::install_github("danwkenn\BELSpatial") # installing the package

# loading the data (Scottish Lip Cancer Data)

# libraries
library(rgdal)          # For readOGR()
library(rgeos)          # For unionSpatialPolygons(), gIntersection(), gBuffer()
library(maptools)       # For unionSpatialPolygons()
library(ggplot2)        # For fortify(), ggplot()
library(readxl)         # For read_excel()
library(magrittr)       # For the pipe operator %>%
library(scales)         # For rescale()
library(dplyr)          # For inner_join(), bind_rows(),                   ### between(), mutate()
library(gridExtra)      # For grid.arrange()
library(tidyr)          # For gather()
library(spdep) # for reading .gal files)
library(MRH) #

map <- readOGR("C:/R dir/Leroux-Empirical-Likelihood-master/Leroux-Empirical-Likelihood-master/scotlip/scotlip.shp", verbose = FALSE)

# Quick plot
plot(map)
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
data1 <- read_excel("C:/R dir/Leroux-Empirical-Likelihood-master/Leroux-Empirical-Likelihood-master/scotlip/scotlip.xlsx") %>% 
  data.frame %>%
  subset(select = -c(1:4, 6:7, 9))
#(data)
#dim(data)


# creating neighbourhood matrix

scot_nb<-read.gal("C:/R dir/Leroux-Empirical-Likelihood-master/Leroux-Empirical-Likelihood-master/scotlip/scotlip.gal", override.id = TRUE)
#class(scot_nb)
W<-nb2mat(scot_nb,style="B")
nblist<-nb2listw(scot_nb)
#creating symmetric neighbourhood matrix for BYM in CARBAYES
rownames(W)<-c()
ind <- upper.tri(W)
W[ind] <- t(W)[ind] 
#isSymmetric(W)
ni<-rowSums(W) # no. of neighbours for each area
R<-diag(ni)
for(i in 1:nrow(R))
{
  R[i,which(W[i,]==1)]<- -1
}


