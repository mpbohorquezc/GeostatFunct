\name{FCOK}
\alias{FCOK}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
FCOK(SFD, newcoords, model, vari = NULL, fill.all = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{SFD}{
%%     ~~Describe \code{SFD} here~~
}
  \item{newcoords}{
%%     ~~Describe \code{newcoords} here~~
}
  \item{model}{
%%     ~~Describe \code{model} here~~
}
  \item{vari}{
%%     ~~Describe \code{vari} here~~
}
  \item{fill.all}{
%%     ~~Describe \code{fill.all} here~~
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
function (SFD, newcoords, model, vari = NULL, fill.all = T) 
{
    if (missing(SFD)) {
        stop("Missing SFD")
    }
    if (missing(newcoords)) {
        stop("Missing new coords")
    }
    if (missing(model)) {
        stop("Missing model")
    }
    if (!inherits(SFD, "SpatFD")) {
        stop("SFD must be an object SpatFD")
    }
    if (!(is.matrix(newcoords) || is.data.frame(newcoords))) {
        stop("Wrong class of newcoords object")
    }
    else if (!all(apply(newcoords, c(1, 2), is.numeric))) {
        stop("Newcoords must be numeric data")
    }
    else if (any(is.na(newcoords))) {
        stop("There is some NA value in newcoords")
    }
    if (missing(vari)) {
        message("Using first variable by default")
    }
    if (missing(fill.all)) {
        message("Using fill.all = TRUE by default")
    }
    if (is.null(vari)) {
        vari = 1
    }
    else if ((is.character(vari) && length(vari) == 1)) {
        if (length(which(names(SFD) == vari)) == 1) {
            vari = which(names(SFD) == vari)
        }
        else if (length(which(names(SFD) == vari)) == 0) {
            stop(paste(vari, "doesn't not exists. Change vari for an existing variable name."))
        }
        else if (length(which(names(SFD) == vari)) == 0) {
            stop("There are more than one variable with the same name")
        }
    }
    if ((is.null(vari) || !(is.numeric(vari) && length(vari) == 
        1))) {
        stop("Wrong class of vari object")
    }
    if (!((isTRUE(fill.all) || isFALSE(fill.all)) && length(fill.all) == 
        1)) {
        stop("Wrong class of fill.all object")
    }
    if (!(inherits(model, "variogramModel") || inherits(model, 
        "list"))) {
        stop("Wrong class of model, model should be of class variogramModel or a list of them (use vgm of gstat package) ")
    }
    else if (inherits(model, "list") && !all(lapply(model, inherits, 
        "variogramModel"))) {
        stop("Wrong class of model, each element of list should be of class variogramModel (use vgm of gstat package)")
    }
    puntaje = list()
    puntajes = list()
    for (k in 1:length(SFD)) {
        puntaje[[k]] = SFD[[k]]$fpca$scores
        rownames(puntaje[[k]]) = SFD[[k]]$coordsname
        puntajes[[k]] = as.data.frame(puntaje[[k]])
        coordinates(puntajes[[k]]) = SFD[[k]]$coords
    }
    colnames(newcoords) = c("x", "y")
    coordinates(newcoords) = ~x + y
    aa = rep(1:length(SFD), lapply(puntajes, ncol))
    bb = unlist(lapply(lapply(puntajes, ncol), seq))
    cc = paste0(c("g=gstat(,paste(colnames(SFD[[", rep("g=gstat(g,paste(colnames(SFD[[", 
        (length(aa) - 1))), aa, rep("]][[\"fpca\"]][[\"harmonics\"]][[\"coefs\"]])[", 
        length(aa)), bb, rep("],", length(aa)), aa, rep(",sep=\".\"),puntajes[[", 
        length(aa)), aa, rep("]][[", length(aa)), bb, rep("]]~1,puntajes[[", 
        length(aa)), aa, rep("]])", length(aa)))
    eval(parse(text = cc))
    g <- gstat(g, model = model, fill.all = fill.all)
    vg <- variogram(g)
    mcl = fit.lmc(vg, g, fit.method = 2)
    plot(vg, model = mcl)
    z = predict(mcl, newdata = newcoords)
    p = 0
    for (i in 1:length(SFD)) {
        p = p + ncol(puntajes[[i]])
    }
    for (i in 1:p) {
        grid.arrange(spplot(z[2 * i - 1], main = "ordinary kriging predictions"), 
            spplot(z[2 * i], main = "ordinary kriging variance"))
    }
    pred = z[1][[1]]
    if (p > 1) {
        for (k in 2:p) {
            pred = cbind(pred, z[2 * k - 1][[1]])
        }
    }
    if (vari == 1) {
        fpred = list()
        for (i in 1:nrow(pred)) {
            fpred[[i]] = SFD[[vari]][["fpca"]][["meanfd"]] + 
                sum((pred[i, 1:ncol(puntajes[[vari]])] * SFD[[vari]][["fpca"]][["harmonics"]]))
        }
    }
    else {
        r = 0
        for (i in 1:(vari - 1)) {
            r = r + ncol(puntajes[[i]])
        }
        fpred = list()
        for (i in 1:nrow(pred)) {
            fpred[[i]] = SFD[[vari]][["fpca"]][["meanfd"]] + 
                sum((pred[i, r:r + ncol(puntajes[[vari]])] * 
                  SFD[[vari]][["fpca"]][["harmonics"]]))
        }
    }
    ret = list(SFD = SFD, model = mcl, fpred = fpred)
    class(ret) = "FCOK"
    return(ret)
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
