\name{scores}
\alias{scores}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Spatial random field of scores
}
\description{
Creates a list of P data-frames, one for each of the P spatial functional variables in the SpatFD-object.  Each data-frame contains one column for each harmonic selected with its respective n_p-score values, and two additional columns for horizontal and vertical spatial coordinates.
}
\usage{
scores(X)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{An object of class SpatFD.
}
}
\details{
For a SpatFD object with P spatial functional variables measured on n_p locations each, and nharm=k_p, p=1,..., P. The scores-function builds a data-frame with n_p rows and k_p+2 columns. The first two columns are the horizontal and vertical spatial coordinates x,y.  The rest of the columns are the score values for each of the first k_p harmonics selected for the p variable, at each of the n_p locations.
}
\value{A list with P data-frames, with dimension n_p \times k_p+2, each. p=1,…,P.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Bohorquez, M., Giraldo, R., & Mateu, J. (2016). Optimal sampling for spatial prediction of functional data. Statistical Methods & Applications, 25(1), 39-54.

Bohorquez, M., Giraldo, R., & Mateu, J. (2016). Multivariate functional random fields: prediction and optimal sampling. Stochastic Environmental Research and Risk Assessment, 31, pages53–70 (2017).
}

\author{Diego Sandoval \email{diasandovalsk@unal.edu.co}.
%%  ~~who you are~~
}
\note{
Functional principal components (FPCA) are applied to the centered data. The scores are scalar second order stationary Random fields, in virtue of the requirements for FPCA. Hence, the covariance function exists, that is, models are bounded, and always have sill and range parameters. Of course,  some models have additional parameters, such as smoothing parameters. From the theoretical perspective, in this case, there is no possibility of non-bounded variogram models. The spatial covariance between two curves is determined by the sum of the covariance between spatial score vectors associated,  see Bohorquez, Giraldo and Mateu 2016 and Bohorquez, Giraldo and Mateu 2021. This covariance can be modeled using the usual packages for geostatistical analysis, such as geoR and gstat.  Finally, the sill of the variogram model for each dimension is bounded for the respective eigenvalue. Actually, an adequate option is to use this eigenvalue as sill and estimate the rest of the parameters.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
PM10=read.table("PM10.txt",head=T,dec=".")
RAMA_Coordenadas=read.table("RAMA_PM10_coordenadas.txt",head=T,row.names=1,dec=",")
estaciones=colnames(PM10)
CoordenadasPM10=RAMA_Coordenadas[estaciones,]
MPM10=as.matrix(PM10,nrow=nrow(PM10),ncol=18,dimnames=c(rownames(PM10[1,4344]),colnames=colnames(PM10)))
PM10spat=SpatFD(MPM10,CoordenadasPM10,basis="Bsplines",nbasis=21,lambda=0.00002,nharm = 2)
NO2=read.table("NO2.txt",head=T,dec=".")
estacionesNO2=read.table("estacionesNO2.txt",head=T,dec=".")
colnames(NO2)=rownames(estacionesNO2)
MNO2=as.matrix(NO2,nrow=4292,ncol=13,dimnames=c(rownames(NO2),colnames=colnames(NO2)))
PM10_MNO2_spat=SpatFD(MNO2, coords = estacionesNO2, basis = "Bsplines", nbasis = 27, lambda = 0.00002, nharm = 4,add=PM10spat)
score_PM10_MNO2=scores(PM10_MNO2_spat)
score_PM10_MNO2[[1]]
score_PM10_MNO2[[2]]
}
## The function is currently defined as
function (X)
{
    if (inherits(X, "SpatFD")) {
        puntaje = list()
        puntajes = list()
        for (j in 1:length(X)) {
            puntaje[[j]] = X[[j]]$fpca$scores
            rownames(puntaje[[j]]) = X[[j]]$coordsnames
            puntajes[[j]] = data.frame(X[[j]]$coords, puntaje[[j]])
        }
    }
    else if (inherits(X, "scores_pred")) {
        puntajes = X$scores_pred
    }
    else {
        stop("Wrong class of X object. It must be of class SpatFD or scores_pred")
    }
    return(puntajes)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{~Functional principal component analysis}
% \keyword{~Spatial process}
% \keyword{~Scores}
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
