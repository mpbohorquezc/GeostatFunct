\name{SpatFD}
\alias{SpatFD}
\title{
Creates univariate and multivariate SpatFD objects.
}
\description{
Creates an object of the class SpatFD from spatial coordinates, and functions or time-series observed at each spatial location.  Time series is a generic term. In fact, observations might be across the frequency or across another spatial dimension such as depth, instead of time.
}
\usage{
SpatFD(data, coords, basis = "Bsplines", nbasis = 4, lambda = 0, nharm = NULL, vp = NULL, name = NULL, add = NULL, ...)
}
\arguments{
  \item{data}{Data must be provided in a data-frame or a matrix where each column corresponds to a location, and the rows are a sequence of data points, that is, the rows are ordered according to time, frequency, depth, …. Data can also be an fd-object from the fda package.
}
  \item{coords}{A data-frame or a matrix with spatial coordinates (x,y). The number of columns in data must coincide with the number of rows in coords for each variable.
}
  \item{basis}{Basis functions. "Fourier" or "Bsplines". By default, "Bsplines".
}
  \item{nbasis}{The number of basis functions.
}
  \item{lambda}{The value of the smoothing parameter.
}
  \item{nharm}{The number of harmonics or eigenfunctions to be reported in the Functional Principal Components results if vp is not given.
}
  \item{vp}{Threshold for the cumulative proportion of explained variance to select the number of components if nharm is not given.
}
\item{name}{A new name for data can be assigned.
}
  \item{add}{
Other variables can be added for spatial multivariate functional prediction, that is, for functional cokriging. It is not necessary that all variables are observed at the same spatial locations.
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{The SpatFD-objects storage the functional data, its parameters, the functional principal component analysis results, and the spatial coordinates for each variable. Each variable has its own functional data, data-frame or matrix and its spatial coordinates file.
}
\value{
For each variable: The functional data and functional principal components linked with spatial coordinates.
}
\references{
Bohorquez, M., Giraldo, R., & Mateu, J. (2016). Optimal sampling for spatial prediction of functional data. Statistical Methods & Applications, 25(1), 39-54.

Bohorquez, M., Giraldo, R., & Mateu, J. (2016). Multivariate functional random fields: prediction and optimal sampling. Stochastic Environmental Research and Risk Assessment, 31, pages53–70 (2017).
}

\author{Diego Sandoval \email{diasandovalsk@unal.edu.co} & Angie Villamil.
}

\note{
1. Although there is no limit for the number of variables for functional cokriging, the real limitation is found on the constraints required to find a valid multivariate covariance model. So, it is highly recommended to apply the parsimony principle.

2. Locations must be in the same region of interest to make sense to include all of them in the same prediction model. However, each variable can be observed in different spatial locations and each can have a different number of observations. There is no limit for the number of variables to be included in this object.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{summary.SpatFD}}
}
\examples{
# Load data
data(BogotaAirQuality)
library(gstat)

# Create an univariate object using 2 nharm
SFD_PM10 <- SpatFD(PM10, coords = coord, basis = "Bsplines", nbasis = 91, lambda = 0.00002, nharm = 2)
View(SFD_PM10)
summary(SFD_PM10)
objects(SFD_PM10)

# Create a bi-variate object using 3 nharm for the first variable and 1 nharm for the second variable
PM10=read.table("PM10.txt",head=T,dec=".")
RAMA_Coordenadas=read.table("RAMA_PM10_coordenadas.txt",head=T,row.names=1,dec=",")
estaciones=colnames(PM10)
CoordenadasPM10=RAMA_Coordenadas[estaciones,]
MPM10=as.matrix(PM10,nrow=nrow(PM10),ncol=18,dimnames=c(rownames(PM10[1,4344]),colnames=colnames(PM10)))
PM10spat=SpatFD(MPM10,CoordenadasPM10,basis="Bsplines",nbasis=21,lambda=0.00002,nharm=3)
summary(PM10spat)
objects(PM10spat)

NO2=read.table("NO2.txt",head=T,dec=".")
estacionesNO2=read.table("estacionesNO2.txt",head=T,dec=".")
nrow(estacionesNO2)
colnames(NO2)=rownames(estacionesNO2)
MNO2=as.matrix(NO2,nrow=4292,ncol=13,dimnames=c(rownames(NO2),colnames=colnames(NO2)))
PM10_MNO2_spatFD=SpatFD(MNO2, coords = estacionesNO2, basis = "Bsplines", nbasis = 27, lambda = 0.00002, nharm = 1,add=PM10spat)
summary(PM10_MNO2_spatFD)
objects(PM10_MNO2_spatFD)
}

\keyword{Functional geostatistics}
\keyword{Spatial data}
