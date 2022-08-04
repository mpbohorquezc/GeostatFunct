\name{SpatFD}
\alias{SpatFD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SpatFD(data, coords, basis = "Bsplines", nbasis = 4, lambda = 0, nharm = NULL, vp = NULL, name = NULL, add = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{coords}{
%%     ~~Describe \code{coords} here~~
}
  \item{basis}{
%%     ~~Describe \code{basis} here~~
}
  \item{nbasis}{
%%     ~~Describe \code{nbasis} here~~
}
  \item{lambda}{
%%     ~~Describe \code{lambda} here~~
}
  \item{nharm}{
%%     ~~Describe \code{nharm} here~~
}
  \item{vp}{
%%     ~~Describe \code{vp} here~~
}
  \item{name}{
%%     ~~Describe \code{name} here~~
}
  \item{add}{
%%     ~~Describe \code{add} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
function (data, coords, basis = "Bsplines", nbasis = 4, lambda = 0, 
    nharm = NULL, vp = NULL, name = NULL, add = NULL, ...) 
{
    if (missing(data)) {
        stop("Missing data")
    }
    if (missing(coords)) {
        stop("Missing coords")
    }
    if (missing(nharm) && missing(vp)) {
        stop("Missing nharm or vp")
    }
    if (!(is.matrix(data) || is.array(data) || is.data.frame(data) || 
        is.fdSmooth(data) || is.fd(data))) {
        stop("Wrong class of data object")
    }
    if (any(is.na(data))) {
        stop("There is some NA value in data")
    }
    if (!(is.matrix(coords) || is.data.frame(coords))) {
        stop("Wrong class of coords object")
    }
    else if (!all(apply(coords, c(1, 2), is.numeric))) {
        stop("Coords must be numeric data")
    }
    else if (any(is.na(coords))) {
        stop("There is some NA value in coords")
    }
    if (is.matrix(data) || is.data.frame(data) || is.array(data)) {
        cx = dim(data)[2]
    }
    else if (is.fdSmooth(data)) {
        cx = dim(data$fd$coefs)[2]
    }
    else if (is.fd(data)) {
        cx = dim(data$coefs)[2]
    }
    fc = dim(coords)[1]
    if (cx != fc) {
        stop("Number of columns of data must be equal to number of rows of coords")
    }
    if (nbasis < nharm) {
        stop("Number of basis must be equal or greater than number of harmonics (nharn)")
    }
    if (!(is.character(basis) && length(basis) == 1)) {
        stop("Wrong class of basis object")
    }
    else if (!(basis == "Fourier" || basis == "Bsplines")) {
        stop("basis not specified")
    }
    if (!(((is.fdSmooth(data) || is.fd(data)) && is.null(nbasis)) || 
        (is.numeric(nbasis) && length(nbasis) == 1))) {
        stop("Wrong class of nbasis object")
    }
    if (!(is.null(nharm) || (is.numeric(nharm) && length(nharm) == 
        1))) {
        stop("Wrong class of nharm object")
    }
    if (!(((is.fdSmooth(data) || is.fd(data)) && is.null(lambda)) || 
        (is.numeric(lambda) && length(lambda) == 1))) {
        stop("Wrong class of lambda object")
    }
    if (!(is.null(vp) || (is.numeric(vp) && length(vp) == 1))) {
        stop("Wrong class of vp object")
    }
    if (!(is.null(add) || inherits(add, "SpatFD"))) {
        stop("Wrong class of add object")
    }
    if (is.null(name)) {
        name = deparse(substitute(data))
    }
    if (!(is.character(name) && length(name) == 1)) {
        stop("Wrong class of name object")
    }
    if (name \%in\% names(add)) {
        stop("Change name, it already exists.")
    }
    if (is.matrix(data) || is.array(data) || is.data.frame(data)) {
        if (missing(basis)) {
            message("Using Bsplines basis by default")
        }
        if (missing(nbasis)) {
            message("Using 4 basis by default")
        }
        if (missing(lambda)) {
            message("Using lambda = 0 by default")
        }
        Mdata = as.matrix(data)
        if (!is.numeric(Mdata)) {
            stop("Object data is not numeric")
        }
        hr <- c(1, nrow(Mdata))
        oplfd <- vec2Lfd(c(1, ncol(Mdata)), hr)
        if (basis == "Bsplines") {
            hourbasis <- create.bspline.basis(hr, nbasis, ...)
        }
        else if (basis == "Fourier") {
            hourbasis <- create.fourier.basis(hr, nbasis, ...)
        }
        data_fdPar <- fdPar(fdobj = hourbasis, Lfdobj = oplfd, 
            lambda)
        data_fdSm <- smooth.basis(argvals = 1:nrow(Mdata), Mdata, 
            data_fdPar)
        data_fd = data_fdSm$fd
        cn = data_fd$fdnames$reps
    }
    else if (is.fdSmooth(data)) {
        data_fdSm = data
        data_fd = data_fdSm$fd
        cn = data_fd$fdnames$reps
    }
    else if (is.fd(data)) {
        data_fd = data
        cn = data_fd$fdnames$reps
    }
    if (!is.null(nharm)) {
        fpca = pca.fd(data_fd, nharm = nharm)
    }
    else if (!is.null(vp)) {
        nh = 1
        repeat {
            fpca = pca.fd(data_fd, nharm = nh)
            if (!sum(fpca$varprop) < vp) {
                break
            }
            nh = nh + 1
        }
    }
    else if (is.null(nharm) && is.null(vp)) {
        fpca = pca.fd(data_fd)
    }
    if (is.null(add)) {
        s = list(list(data = data, coords = coords, coordsnames = cn, 
            data_fd = data_fd, fpca = fpca, variable_name = name))
        names(s) = name
        class(s) = "SpatFD"
    }
    else if (inherits(add, "SpatFD")) {
        s = list(list(data = data, coords = coords, coordsnames = cn, 
            fpca = fpca, variable_name = name))
        names(s) = name
        s = append(add, s)
        class(s) = "SpatFD"
    }
    return(s)
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
