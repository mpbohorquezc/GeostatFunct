\name{recons_fd}
\alias{recons_fd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
This is an internal function for functional kriging and cokriging. To perform the linear combinations to obtain functional kriging and cokriging.
}
Once optimization process is finished and scores or lambda are obtained, this function builds the prediction, performing the linear combination between coefficients and basis functions.
}
\usage{
recons_fd(X)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%~~A KS_pred object obtained with function KS_scores_lambdas~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Spatial functional predictions. The fd object based on the two functional kriging methods described in Bohorquez, M., Giraldo, R., & Mateu, J. (2016).
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~
Bohorquez, M., Giraldo, R., & Mateu, J. (2016). Optimal sampling for spatial prediction of functional data. Statistical Methods & Applications, 25(1), 39-54.

Bohorquez, M., Giraldo, R., & Mateu, J. (2016). Multivariate functional random fields: prediction and optimal sampling. Stochastic Environmental Research and Risk Assessment, 31, pages53–70 (2017).

Bohorquez M.; Giraldo R. and Mateu J. Spatial prediction and optimal sampling of functional data in Geostatistical Functional Data Analysis: Theory and Methods (2021).  John Wiley  Sons, Chichester, UK. ISBN: 978-1-119-38784-8. https://www.wiley.com/en-us/Geostatistical+Functional+Data+Analysis-p-9781119387848. ~
}
\author{
Diego Sandoval \email{diasandovalsk@unal.edu.co}.
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
curves_PM10_l <- recons_fd(KS_SFD_PM10_l)
plot(curves_PM10_l)
#method = "scores"
curves_PM10_sc <- recons_fd(KS_SFD_PM10_sc)
plot(curves_PM10_sc)
}

\keyword{Spatial prediction}
\keyword{Eigenfunctions}
