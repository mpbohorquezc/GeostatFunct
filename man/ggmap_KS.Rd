\name{ggmap_KS}
\alias{ggmap_KS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ggmap_KS(KS, map_path, window_time = NULL, method = "lambda", map_n = 5000, zmin = NULL, zmax = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{KS}{
%%     ~~Describe \code{KS} here~~
}
  \item{map_path}{
%%     ~~Describe \code{map_path} here~~
}
  \item{window_time}{
%%     ~~Describe \code{window_time} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
}
  \item{map_n}{
%%     ~~Describe \code{map_n} here~~
}
  \item{zmin}{
%%     ~~Describe \code{zmin} here~~
}
  \item{zmax}{
%%     ~~Describe \code{zmax} here~~
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
function (KS, map_path, window_time = NULL, method = "lambda", 
    map_n = 5000, zmin = NULL, zmax = NULL) 
{
    map <- readOGR(map_path)
    newcoords <- sp::spsample(map, map_n, type = "regular")
    newcoords <- as.data.frame(newcoords)
    colnames(newcoords) <- c("X", "Y")
    KS_SFD <- KS_scores_lambdas(KS$SFD, newcoords, model = KS$model, 
        method = method, name = KS$name)
    SFD <- recons_fd(KS_SFD)
    if (is.null(window_time)) {
        times <- SFD$basis$rangeval[1]
    }
    else if (!(all(window_time >= SFD$basis$rangeval[1]) && all(window_time <= 
        SFD$basis$rangeval[2]))) {
        stop(paste("window_time is out of bounds: Must be some value(s) between ", 
            SFD$basis$rangeval[1], "and ", SFD$basis$rangeval[2]))
    }
    else {
        times <- sort(window_time)
    }
    eval <- eval.fd(times, SFD)
    melt_s <- suppressWarnings(melt(eval))
    melt_s$X2 <- as.factor(melt_s$X2)
    melt_s$X <- as.factor(melt_s$X2)
    levels(melt_s$X) <- newcoords$X
    melt_s$Y <- as.factor(melt_s$X2)
    levels(melt_s$Y) <- newcoords$Y
    names(melt_s) = c("Time", "Prediction", "Value", "X", "Y")
    melt_s$Time <- as.factor(melt_s$Time)
    graf <- list()
    if (is.null(zmin)) {
        zmin = min(melt_s$Value)
    }
    if (is.null(zmax)) {
        zmax = max(melt_s$Value)
    }
    for (i in 1:nlevels(melt_s$Time)) {
        melt_s_2 <- melt_s[melt_s$Time == i, ]
        graf[[i]] <- plot_ly(x = as.numeric(as.character(melt_s_2$X)), 
            y = as.numeric(as.character(melt_s_2$Y)), z = melt_s_2$Value, 
            type = "heatmap", colorbar = list(title = "Prediction"), 
            reversescale = T, zmin = zmin, zmax = zmax) \%>\% layout(title = paste("Prediction - Time = ", 
            times[i]), xaxis = list(showticklabels = FALSE), 
            yaxis = list(showticklabels = FALSE), scene = list(aspectration = list(x = 1, 
                y = 1)))
    }
    return(graf)
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
