\name{FD_optimal_design}
\alias{FD_optimal_design}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Optimal Spatial Design For Functional Data
}
\description{
Given a variogram model and a set of points in which we want to predict certain variable optimally, this function finds where must be placed the stations in which the information will be collected.
}
\usage{
FD_optimal_design(k, s0, model, fixed_stations = NULL,
                  method = "lambda", grid = NULL,
                  map = NULL, nharm = NULL, plt = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{k}{The number of new stations that are going to be located.}
  \item{s0}{matrix, array, data.frame, or SpatialPoints object that contains the coordinates of the points in which we want to optimally predict using functional kriging. Each row corresponds to each point.
  If you want to predict optimally in a certain area, you must create a random sample large enought of points of that area and pass it to this function.}
  \item{model}{A VariogramModel from gstat package or a list of VariogramModels if a distinct model is going to be used for each harmonic.}
  \item{fixed_stations}{If there are already some stations on the field that are not goint to be removed, here should be passed their coordinates. The object must be of class matrix, array, data.frame, SpatialPoints,SpatFD, or NULL if there are not fidex stations.}
  \item{method}{Functional kriging method that will be used. Currently available "lambda" and "scores". See details bellow.}
  \item{grid}{Coordinates in which the new stations can be located. grid must be of type matrix, array, data.frame, SpatialPoints. If you don't pass grid, you must pass map; if you pass both, map will be used.}
  \item{map}{Spatial object from sp package such as Line, Lines, Polygon, SpatialPolygons, SpatialGrid or SpatialPixels in which the new stations will be located. This object will also be used for creating the plot.}
  \item{nharm}{Number of harmonics of the functional principal components that will be used in the prediction.}
  \item{plt}{TRUE or FALSE. If TRUE, a nice ggplot2 plot with the result will be generated. Se example bellow.}
}
\details{

"scores" method consists of...

On the other hand "lambda" method...

See Bohorquez, M., Giraldo, R., & Mateu, J. (2016).

}
\value{

The function returns an OptimalSpatialDesign object that is a list with the following elements:

  \item{new_stations}{matrix array object with the coordinates of the new stations.}
  \item{fixed_stations}{matrix array object with the coordinates of the fixed stations.}
  \item{plot}{ggplot2 plot.}
  
}
\references{
Bohorquez, M., Giraldo, R., & Mateu, J. (2016). Optimal sampling for spatial prediction of functional data. Statistical Methods & Applications, 25(1), 39-54.
}
\author{
Nathaly Vergel Serrano and Samuel Sánchez Gutiérrez.
National University of Colombia, Statistics Department.
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~

}
\examples{
s0 <- cbind(2*runif(100),runif(100)) # random coordinates on (0,2)x(0,1)
fixed_stations <- cbind(2*runif(4),runif(4))
x_grid <- seq(0,2,length = 30)
y_grid <- seq(0,1,length = 30)
grid <- cbind(rep(x_grid,each = 30),rep(y_grid,30))
model  <- gstat::vgm(psill = 5.665312,
                  model = "Exc",
                  range = 8000,
                  kappa = 1.62,
                  add.to = vgm(psill = 0.893,
                               model = "Nug",
                               range = 0,
                               kappa = 0))
#
FD_optimal_design(k = 10, s0 = s0, model = model,
                  grid = grid, nharm = 2, plt = T,
                  fixed_stations = fixed_stations) -> OSD
#
OSD
class(OSD)
plot(OSD)
}

