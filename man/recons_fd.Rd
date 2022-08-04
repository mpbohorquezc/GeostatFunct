\name{recons_fd}
\alias{recons_fd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
recons_fd(X)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{X} here~~
}
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
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
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
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
