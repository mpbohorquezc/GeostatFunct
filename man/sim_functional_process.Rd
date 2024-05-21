\name{sim_functional_process}
\alias{sim_functional_process}
\title{
Simulation of unconditional or conditional functional spatial process.
}
\description{
Given a variogram model, this functions simulates several realizations of a functional spatial process. This simulation can be conditioned to observed data.
}
\usage{
sim_functional_process <- function(nsims,variograms,
                                   nbasis,coords,
                                   data = NULL,
                                   data_coords = NULL,
                                   basis = NULL,
                                   mu = NULL,
                                   L = NULL)
}
\arguments{
  \item{nsims}{Integer giving the number of curves to simulate.}
  \item{variograms}{\code{gstat::gstatVariogram} or \code{list} of them giving the variogram model for each score. If only one is provided, it will be recycled.}
  \item{nbasis}{Integer giving the number of basis functions on which the process is going to be projected.}
  \item{coords}{Gridded \code{sp::SpatialPoints} or \code{sp::SpatialPixels}, or array coordinates of the curves that are going to be simulated.}
  \item{data}{\code{fda::fd} object containing the observed curves for conditional simulation. If data is not provided, inconditional simulation is performed.}
  \item{data_coords}{\code{sp::SpatialPoints} or array coordinates of the observed data.}
  \item{basis}{Character giving the basis of functions (only for inconditional simulation) (nbasis must be provided).}
  \item{mu}{\code{fda::fd} object of the mean function of the process, default is zero. Only used in unconditional simulation.}
  \item{L}{Limits of the symetric interval centered on zero that is the domain of the basis that is going to be created in unconditional simulation case.}
}
\details{
When \code{data} is passed, conditional simulation is performed. That means that each simulated realization of the process interpolated the observed curves in \code{data}. If \code{data} is \code{NULL}, the realizations of the process are simulated without imterpolation restrictions.
}
\value{
A list of \code{nsims} \code{fd} objects each one with as much curves as points are in \code{coords}.
}
\references{
}
\author{
Samuel Sánchez Gutiérrez \email{ssanchezgu@unal.edu.co}.
}

\section{Warning}{
}

\seealso{\code{\link{generate_basis}}
}
\examples{
\donttest{
}
}
