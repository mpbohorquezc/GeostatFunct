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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (X)
{
    if (inherits(X, "KS_pred")) {
        vari = X$name
        if (length(X) == 4) {
            if (inherits(X[[2]], "scores_pred")) {
                a = scores(X[[2]])
                mean_coef = X$SFD[[vari]]$fpca$meanfd$coefs
                nr = nrow(a)
                mean_coef = matrix(rep(mean_coef, nr), ncol = nr)
                coef_scores = (X$SFD[[vari]]$fpca$harmonics$coef \%*\%
                  t(a)) + mean_coef
                result_scores = fd(coef_scores, X$SFD[[vari]]$fpca$harmonics$basis)
                result = result_scores
            }
            else {
                b = as.matrix(X[[2]]$lambda_pred)
                mean_coef = X$SFD[[vari]]$fpca$meanfd$coefs
                nr = ncol(b)
                mean_coef = matrix(rep(mean_coef, nr), ncol = nr)
                coef_lambda = (X$SFD[[vari]]$data_fd$coefs \%*\%
                  b) + mean_coef
                result_lambda = fd(coef_lambda, X$SFD[[vari]]$data_fd$basis)
                result = result_lambda
            }
        }
        else {
            a = scores(X[[2]])
            mean_coef = X$SFD[[vari]]$fpca$meanfd$coefs
            nr = nrow(a)
            mean_coef = matrix(rep(mean_coef, nr), ncol = nr)
            coef_scores = (X$SFD[[vari]]$fpca$harmonics$coef \%*\%
                t(a)) + mean_coef
            result_scores = fd(coef_scores, X$SFD[[vari]]$fpca$harmonics$basis)
            b = as.matrix(X[[3]]$lambda_pred)
            mean_coef = X$SFD[[vari]]$fpca$meanfd$coefs
            nr = nrow(a)
            mean_coef = matrix(rep(mean_coef, nr), ncol = nr)
            coef_lambda = (X$SFD[[vari]]$data_fd$coefs \%*\% b) +
                mean_coef
            result_lambda = fd(coef_lambda, X$SFD[[vari]]$data_fd$basis)
            result = list(fd_scores = result_scores, fd_lambda = result_lambda)
        }
    }
    else {
        stop("Wrong class KS_pred")
    }
    return(result)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{~Spatial prediction}
% \keyword{~Eigenfunctions}
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
Footer
© 2023 GitHub, Inc.
Footer navigation
Terms
Privacy
Security
Status
Docs
Contact GitHub
Pricing
API
Training
Blog
About

