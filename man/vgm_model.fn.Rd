\name{vgm_model.fn}
\alias{vgm_model.fn}
\title{
Functional kriging Variance
}
\description{
This is an internal function of the package that is used by \code{FD_optimal_design} function to calculate a value directly proportional to the variance of functional kriging.
}

\usage{

vgm_model.fn(points, s0, model,
            method_i = "lambda", grid,
            fixed_stations)
}
\arguments{
  \item{points}{current position of the new stations.}
  \item{s0}{\code{matrix}, \code{array}, \code{data.frame}, or \code{SpatialPoints} object that contains the coordinates of the points in which we want to optimally predict using functional kriging. Each row corresponds to each point.}
  \item{model}{A \code{VariogramModel} from \code{gstat} package or a list of \code{VariogramModels} objects if a distinct model is going to be used for each harmonic.}
  \item{method_i}{Functional kriging method that will be used. Currently available "lambda" and "scores". See details bellow.}
  \item{grid}{Coordinates in which the new stations can be located. \code{grid} must be of type \code{matrix}, \code{array}, \code{data.frame}, \code{SpatialPoints}.}
  \item{fixed_stations}{If there are already some stations on the field that are not goint to be removed, here should be passed their coordinates. The object must be of class \code{matrix}, \code{array}, \code{data.frame}, \code{SpatialPoints}, \code{SpatFD}, or \code{NULL} if there are no fidex stations.}

}
%\details{

%}
\value{
A scalar. For 'scores' method, the simple kriging variance, for 'lambda' method a value proportional to the kriging variance.
}
\references{
Bohorquez, M., Giraldo, R., & Mateu, J. (2016). Optimal sampling for spatial prediction of functional data. Statistical Methods & Applications, 25(1), 39-54.
}
\author{
Nathaly Vergel Serrano \email{nvergel@unal.edu.co} & Samuel Sánchez Gutiérrez \email{ssanchezgu@unal.edu.co}.
}

\section{warning}{
  Please do not use this function directly.
}

\seealso{
\code{\link[SpatFD]{FD_scores_lambda}}
}
\examples{
# See FD_optimal_design examples
}

