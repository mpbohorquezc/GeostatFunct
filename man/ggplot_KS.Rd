\name{ggplot_KS}
\alias{ggplot_KS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ggplot_KS(KS, show.varpred = T, main = "Functional Data", main2 = "Functional Data", ylab = "Value", xlab = "Time", ndigits = 3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{KS}{
%%     ~~Describe \code{KS} here~~
}
  \item{show.varpred}{
%%     ~~Describe \code{show.varpred} here~~
}
  \item{main}{
%%     ~~Describe \code{main} here~~
}
  \item{main2}{
%%     ~~Describe \code{main2} here~~
}
  \item{ylab}{
%%     ~~Describe \code{ylab} here~~
}
  \item{xlab}{
%%     ~~Describe \code{xlab} here~~
}
  \item{ndigits}{
%%     ~~Describe \code{ndigits} here~~
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
function (KS, show.varpred = T, main = "Functional Data", main2 = "Functional Data", 
    ylab = "Value", xlab = "Time", ndigits = 3) 
{
    if (missing(KS)) {
        stop("Missing KS")
    }
    if (!inherits(KS, "KS_pred")) {
        stop("KS must be an object KS_pred")
    }
    SFD <- recons_fd(KS)
    if (length(SFD) != 2) {
        times <- SFD$basis$rangeval[1]:SFD$basis$rangeval[2]
        eval <- eval.fd(times, SFD)
        melt_s <- suppressWarnings(melt(eval))
        melt_s$X2 <- as.factor(melt_s$X2)
        if (show.varpred) {
            if (inherits(KS[[2]], "scores_pred")) {
                levels(melt_s$X2) <- paste(levels(melt_s$X2), 
                  " \n(var_pred = ", round(KS[[2]]$scores_varpred$VTotal, 
                    ndigits), ")", sep = "")
                if (missing(main)) {
                  main <- "Functional Data - Scores Method"
                }
            }
            else {
                levels(melt_s$X2) <- paste(levels(melt_s$X2), 
                  " \n(var_pred = ", round(KS[[2]]$lambda_varpred$VTotal, 
                    ndigits), ")", sep = "")
                if (missing(main)) {
                  main <- "Functional Data - Lambda Method"
                }
            }
        }
        names(melt_s) = c("Time", "Prediction", "Value")
        graf = ggplot(melt_s, aes(x = Time, y = Value, col = Prediction)) + 
            geom_line() + labs(title = main) + labs(x = xlab, 
            y = ylab) + theme_minimal()
        return(graf)
    }
    else {
        graf <- list()
        mainl <- list(main, main2)
        if (missing(main)) {
            mainl[[1]] <- "Functional Data - Scores Method"
        }
        if (missing(main2)) {
            mainl[[2]] <- "Functional Data - Lambda Method"
        }
        for (i in 1:2) {
            times <- SFD[[i]]$basis$rangeval[1]:SFD[[i]]$basis$rangeval[2]
            eval <- eval.fd(times, SFD[[i]])
            melt_s <- suppressWarnings(melt(eval))
            melt_s$X2 <- as.factor(melt_s$X2)
            if (show.varpred) {
                if (inherits(KS[[1 + i]], "scores_pred")) {
                  levels(melt_s$X2) <- paste(levels(melt_s$X2), 
                    " \n(var_pred = ", round(KS[[1 + i]]$scores_varpred$VTotal, 
                      ndigits), ")", sep = "")
                }
                else {
                  levels(melt_s$X2) <- paste(levels(melt_s$X2), 
                    " \n(var_pred = ", round(KS[[1 + i]]$lambda_varpred$VTotal, 
                      ndigits), ")", sep = "")
                }
            }
            names(melt_s) = c("Time", "Prediction", "Value")
            graf[[i]] = ggplot(melt_s, aes(x = Time, y = Value, 
                col = Prediction)) + geom_line() + labs(title = mainl[[i]]) + 
                labs(x = xlab, y = ylab) + theme_minimal()
        }
        return(graf)
    }
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
