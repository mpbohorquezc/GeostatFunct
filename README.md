SpatFD
======================
Functional geostatistics: univariate and multivariate functional spatial prediction

![version](https://img.shields.io/badge/version-0.1.0-blue)
[![Licence](https://img.shields.io/badge/licence-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)

The R package *SpatFD* carries out functional kriging  and optimal sampling for spatial prediction of functional data. The methods are based on XAs and YAs  that extend the framework of  spatial prediction and optimal sampling from scalar to functional data. SpatFD is based on the KHH expansion that allows to represent the observed functions in terms of its empirical functional principal components. Based on this approach, the functional auto-covariances and cross-covariances required for  spatial functional predictions and optimal sampling, are completely determined by the sum of the spatial auto-covariances and cross-covariances of the respective score components.  The package provides new classes of data and  functions for modeling spatial dependence structure among curves. The spatial prediction of curves at unsampled locations can be carried out using two types of predictors, and both of them report, the respective variances of the prediction error.  In addition, there is a function for the determination of spatial locations sampling configuration that ensures minimum variance of spatial functional prediction. There are also two functions for plotting predicted curves at each location and mapping the surface at each time point, respectively.

## Installation
You can install the **development** version from [Github](https://github.com/mpbohorquezc/GeostatFunct).
```s
install.packages("devtools")
devtools::install_github("mpbohorquezc/GeostatFunct")
```

## Overview
The objects class and usage in the different functions are listed.

- `SpatFD` return an object class 'SpatFD' which is used in functional kriging and to obtain the spatial random field of scores.
- `KS_scores_lambdas` return an object class 'KS_pred' to use in the linear combinations to obtain functional kriging and plots with `ggplot_KS` and `ggmap_KS`.
- `FD_optimal_design` return an object class 'OptimalSpatialDesign' that can be used with print.

## Example of use
```s
library(SpatFD)
library(gstat)

# load data and coordinates
data(AirQualityBogota)

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
class(KS_SFD_PM10_l)

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
         map_path = map,
         window_time = c(3500),
         zmin = 25,
         zmax = 100)

ggmap_KS(KS_SFD_PM10_both,
         map_path = map,
         window_time = c(2108),
         method = "lambda",
         zmin = 50,
         zmax = 120)

ggmap_KS(KS_SFD_PM10_both,
         map_path = map,
         window_time = c(5108,5109,5110),
         method = "scores",
         zmin = 50,
         zmax = 120)
```

## For more information
You can read:

  - https://onlinelibrary.wiley.com/doi/book/10.1002/9781119387916
  - https://link.springer.com/article/10.1007/s10260-015-0340-9
  - https://link.springer.com/article/10.1007/s00477-016-1266-y

## License
This package is free and open source software, licensed under GPL-3.

## References
 * Bohorquez, M., Giraldo, R., & Mateu, J. (2016). Optimal sampling for spatial prediction of functional data. Statistical Methods & Applications, 25(1), 39-54.
 * Bohorquez, M., Giraldo, R., & Mateu, J. (2016). Multivariate functional random fields: prediction and optimal sampling. Stochastic Environmental Research and Risk Assessment, 31, pages53–70 (2017).
 * Bohorquez M.; Giraldo R. and Mateu J. Spatial prediction and optimal sampling of functional data in Geostatistical Functional Data Analysis: Theory and Methods (2021).  John Wiley  Sons, Chichester, UK. ISBN: 978-1-119-38784-8. https://www.wiley.com/en-us/Geostatistical+Functional+Data+Analysis-p-9781119387848.


