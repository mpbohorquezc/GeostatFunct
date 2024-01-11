\name{COKS_scores_lambdas}
\alias{COKS_scores_lambdas}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Functional CoKriging
}
\description{
Linear Spatial functional prediction. Two predictors are possible: scores or lambda.
}
\usage{
COKS_scores_lambdas(SFD, newcoords, model, method = "scores", fill.all=T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{SFD}{
object of class 'SpatFD'.
}
  \item{newcoords}{
The N × 2 matrix or data.frame with the spatial coordinates corresponding to the N prediction locations.
}
  \item{model}{
The scores variogram model. A variogramModel object. See gstat package.
}  
\item{method}{
Prediction method: "scores"
}
  \item{fill.all}{
gstat function parameter. If there are more than 1 score vector and not all models or a valid and complete linear model of coregionalization is given, fill all of the direct and cross variogram model with the only model given.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a 'COKS_pred' object with functional kriging
}
\references{
Bohorquez, M., Giraldo, R., & Mateu, J. (2016). Optimal sampling for spatial prediction of functional data. Statistical Methods & Applications, 25(1), 39-54.

Bohorquez, M., Giraldo, R., & Mateu, J. (2016). Multivariate functional random fields: prediction and optimal sampling. Stochastic Environmental Research and Risk Assessment, 31, pages53–70 (2017).

Bohorquez M.; Giraldo R. and Mateu J. Spatial prediction and optimal sampling of functional data in Geostatistical Functional Data Analysis: Theory and Methods (2021).  John Wiley  Sons, Chichester, UK. ISBN: 978-1-119-38784-8. https://www.wiley.com/en-us/Geostatistical+Functional+Data+Analysis-p-9781119387848.
}
\author{
Valeria Bejarano <vbejaranos@unal.edu.co>
}
\seealso{
\code{\link{SpatFD}},\code{\link{summary.COKS_pred}}
}
\examples{
data(COKMexico)
SFD_PM10_NO2 <- SpatFD(Mex_PM10, coords = coord_PM10, basis = "Fourier", 
nbasis = 21, lambda = 0.000001, nharm = 2)
SFD_PM10_NO2 <- SpatFD(NO2, coords = coord_NO2, basis = "Fourier", 
nbasis = 27, lambda = 0.000001, nharm = 2,
                      add = SFD_PM10_NO2)
model1 <- gstat::vgm(647677.1,"Gau",23317.05)
model1 <- gstat::vgm(127633,"Wav",9408.63, add.to = model1)
newcoords <- data.frame(x = 509926,y = 2179149)
COKS_scores_lambdas(SFD_PM10_NO2,newcoords,model1)
}
\keyword{Cokriging}
\keyword{Spatial functional prediction}