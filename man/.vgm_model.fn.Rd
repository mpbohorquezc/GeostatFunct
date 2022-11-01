\name{.vgm_model.fn}
\alias{.vgm_model.fn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Functional kriging Variance
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This is an internal function of the package that is used by \link{OptimalSpatialDesign} function to calculate a value directly proportional to the variance of functional kriging.
}

\usage{

.vgm_model.fn(points, s0, model,
            method_i = "lambda", grid, 
            fixed_stations)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{points}{current position of the new stations.}
  \item{s0}{\code{matrix}, \code{array}, \code{data.frame}, or \code{SpatialPoints} object that contains the coordinates of the points in which we want to optimally predict using functional kriging. Each row corresponds to each point.}
  \item{model}{A \code{VariogramModel} from \code{gstat} package or a list of \code{VariogramModels} objects if a distinct model is going to be used for each harmonic.}
  \item{method_i}{Functional kriging method that will be used. Currently available "lambda" and "scores". See details bellow.}
  \item{grid}{Coordinates in which the new stations can be located. \code{grid} must be of type \code{matrix}, \code{array}, \code{data.frame}, \code{SpatialPoints}.}
  \item{fixed_stations}{If there are already some stations on the field that are not goint to be removed, here should be passed their coordinates. The object must be of class \code{matrix}, \code{array}, \code{data.frame}, \code{SpatialPoints}, \code{SpatFD}, or \code{NULL} if there are no fidex stations.}
  %\item{x}{
%%     ~~Describe \code{x} here~~
%}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Nathaly Vergel Serrano and Samuel Sánchez Gutiérrez.
National University of Colombia, Statistics Department.
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\section{warning}{
  Please do not use this function directly.
}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
