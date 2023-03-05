rm(list=ls())

library(gstat)
library(sp)
library(fda)
-library(reshape)
-library(geoR)
library(ggplot2)
library(rgdal)
library(sf)
library(plotly)

source("SpatFD.R")
source("outputs.R")
source("scores.R")
source("recons_fd.R")
source("KS_scores_lambdas.R")
source("ggplot_KS.R")
source("ggmap_KS.R")
source("OSD_scores_lambda.R")

#data(AirQualityBogota)
# load data and coordinates
PM10 = read.table("PartMat10.txt",head=T,dec=",")
coord = read.table("Coords_DAMABog.txt",dec=",",sep="\t",header=T)

#s_0 nonsampled location. It could be data.frame or matrix and one or more locations of interest
newcoorden=data.frame(X=seq(93000,105000,len=100),Y=seq(97000,112000,len=100))
#newcoorden=data.frame(X=110000,Y=126000)
#newcoorden=matrix(c(110000.23,109000,109500,130000.81,129000,131000),nrow=3,ncol=2,byrow=T)

# Building the SpatFD object
SFD_PM10 <- SpatFD(PM10, coords = coord[, -1], basis = "Bsplines", nbasis = 17,norder=5, lambda = 0.00002, nharm=3)
summary(SFD_PM10)

#Semivariogram models for each spatial random field of scores
modelos <- list(vgm(psill = 2634000, "Exp", range = 2103.25, nugget =  0),
                vgm(psill = 101494.96, "Exp", range = 1484.57, nugget = 0),
                vgm(psill =53673, "Exp", range = 42406, nugget =  0))

#Functional kriging. Functional spatial prediction at each location of interest
#method = "lambda"
#Computation of lambda_i
KS_SFD_PM10_l <- KS_scores_lambdas(SFD_PM10, newcoorden ,method = "lambda", model = modelos)
#method = "scores"
#Simple kriging of scores
KS_SFD_PM10_sc <- KS_scores_lambdas(SFD_PM10, newcoorden, method = "scores", model = modelos)
#method = "both"
KS_SFD_PM10_both <- KS_scores_lambdas(SFD_PM10, newcoorden, method = "both", model = modelos)

summary(KS_SFD_PM10_l)
summary(KS_SFD_PM10_sc)
summary(KS_SFD_PM10_both)

#Linear combinations among weigths predictions and eigenfunctions
recons_fd(KS_SFD_PM10_l)
recons_fd(KS_SFD_PM10_sc)
recons_fd(KS_SFD_PM10_both)

#Curve and variance prediction plots
ggplot_KS(KS_SFD_PM10_l)
ggplot_KS(KS_SFD_PM10_l, show.varpred = F) 
ggplot_KS(KS_SFD_PM10_sc)
ggplot_KS(KS_SFD_PM10_sc, show.varpred = F) 
#Curve and variance prediction for both methods
PlotKS=ggplot_KS(KS_SFD_PM10_both,
          main = "Plot 1 - Using Scores",
          main2 = "Plot 2 - Using Lambda",
          ylab = "PM10")
PlotKS[[1]]
PlotKS[[2]]

#Smoothed prediction maps for the given specific times 
ggmap_KS(KS_SFD_PM10_l,
         map_path = "Bogota.shp",
         window_time = c(3500),
         zmin = 25,
         zmax = 100)

ggmap_KS(KS_SFD_PM10_both,
         map_path = "Bogota.shp",
         window_time = c(2108),
         method = "lambda",
         zmin = 50,
         zmax = 120)

ggmap_KS(KS_SFD_PM10_both,
         map_path = "Bogota.shp",
         window_time = c(5108,5109,5110),
         method = "scores",
         zmin = 50,
         zmax = 120)




