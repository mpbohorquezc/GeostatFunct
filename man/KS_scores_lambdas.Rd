\name{KS_scores_lambdas}
\alias{KS_scores_lambdas}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Functional Kriging
}
\description{
Linear Spatial functional prediction. Two predictors are possible: scores or lambda.
}
\usage{
KS_scores_lambdas(SFD, newcoords, model, method = "lambda", name = NULL, fill.all = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{SFD}{
SpatFD object.
}
  \item{newcoords}{
The N × 2 matrix or data-frame with the spatial coordinates corresponding to prediction locations.
}
  \item{model}{
The scores variogram model. A for variogramModel object. See gstat package.
}
  \item{method}{
Prediction method: "lambda" or "scores". By default "lambda". "lambda" option corresponds to the predictor $\breve{\chi}_{\bm s_0}(t)$, given by
\[\breve{\chi}_{\bm s_0}(t)=\sum\limits_{i=1}^{n}\lambda_i\chi_{\bm s_i}(t)\]
and weigths are found such that minimize
$\left\|\bm\chi_{\bm s_0}(t)-\breve{\chi}_{\bm s_0}(t)\right\|^2$

"scores" method performs kriging or cokriging directly on the scores and predictions are used to build the functional prediction
}
  \item{name}{
The variable to predict in SpatFD object. By default, the predictions is performed for the first variable in the SpatFD object.
}
  \item{fill.all}{
gstat function parameter. If there are more than 1 score vector and not all models or a valid and complete linear model of coregionalization is given, fill all of the direct and cross variogram model with the only model given.
}
}
Bohorquez M.; Giraldo R. and Mateu J. Spatial prediction and optimal sampling of functional data in Geostatistical Functional Data Analysis: Theory and Methods (2021).
}
\value{
%%  ~Functional kriging weights:
%%  weights (lambda) using the first method
%%  and kriging score predictions using the second method in Bohorquez, M., Giraldo, R., & Mateu, J. (2016).
%%  Variance for each spatial functional prediction
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Bohorquez, M., Giraldo, R., & Mateu, J. (2016). Optimal sampling for spatial prediction of functional data. Statistical Methods & Applications, 25(1), 39-54.

Bohorquez, M., Giraldo, R., & Mateu, J. (2016). Multivariate functional random fields: prediction and optimal sampling. Stochastic Environmental Research and Risk Assessment, 31, pages53–70 (2017).

Bohorquez M.; Giraldo R. and Mateu J. Spatial prediction and optimal sampling of functional data in Geostatistical Functional Data Analysis: Theory and Methods (2021).  John Wiley  Sons, Chichester, UK. ISBN: 978-1-119-38784-8. https://www.wiley.com/en-us/Geostatistical+Functional+Data+Analysis-p-9781119387848.
}
\author{
Diego Sandoval \email{diasandovalsk@unal.edu.co} & Angie Villamil.
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
PM10 = read.table("PartMat10.txt",head=T,dec=",")
coord = read.table("Coords_DAMABog.txt",dec=",",sep="\t",header=T,row.names=1)
#newcoorden=matrix(c(110000.23,125000,109500,131000.81,108000,134000),nrow=3,ncol=2,byrow=T)
newcoorden=data.frame(X=110000,Y=125000)
# Recibir los datos, suavizarlos y ACP
SFD_PM10 <- SpatFD(PM10, coords = coord, basis = "Bsplines", nbasis = 17,norder=5, lambda = 0.00002, nharm=3)

#Variogram model for each component
modelos <- list(vgm(psill = 2634000, "Exp", range = 2103.25, nugget =  0),
                vgm(psill = 101494.96, "Exp", range = 1484.57, nugget = 0),
                vgm(psill =53673, "Exp", range = 42406, nugget =  0))

#Genera los scores y los lambdas para predecir en nuevas coordenadas

#method = "lambda"
KS_SFD_PM10_l <- KS_scores_lambdas(SFD_PM10, newcoorden ,method = "lambda", model = modelos)
#method = "scores"
KS_SFD_PM10_sc <- KS_scores_lambdas(SFD_PM10, newcoorden, method = "scores", model = modelos)
#method = "both"
KS_SFD_PM10_both <- KS_scores_lambdas(SFD_PM10, newcoorden, method = "both", model = modelos)

# Retorna a SFD_PM10, scores, lambdas y las varianzas de predicción

#method = "lambda"
KS_SFD_PM10_l$KS_lambda$lambda_pred
KS_SFD_PM10_l$KS_lambda$lambda_varpred

#method = "scores"
KS_SFD_PM10_sc$KS_scores$scores_pred
KS_SFD_PM10_sc$KS_scores$scores_varpred

#method = "both"
KS_SFD_PM10_both$KS_lambda$lambda_pred
KS_SFD_PM10_both$KS_lambda$lambda_varpred
KS_SFD_PM10_both$KS_scores$scores_pred
KS_SFD_PM10_both$KS_scores$scores_varpred

summary(KS_SFD_PM10_both)

PlotKS=ggplot_KS(KS_SFD_PM10_both,
          main = "Plot 1 - Using Scores",
          main2 = "Plot 2 - Using Lambda",
          ylab = "PM10")
PlotKS[[1]]
PlotKS[[2]]

ggplot_KS(KS_SFD_PM10_l, show.varpred = F)

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (SFD, newcoords, model, method = "lambda", name = NULL,
    fill.all = NULL)
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
    if (missing(name)) {
        message("Using first variable by default")
    }
    if (missing(fill.all)) {
        message("Using fill.all = TRUE by default")
    }
    if (missing(method)) {
        message("Using method = 'lambda' by default")
    }
    if (!(is.character(method) && length(method) == 1)) {
        stop("Wrong class of method object")
    }
    else if (!(method == "lambda" || method == "scores" || method ==
        "both")) {
        stop("method must be one of the following: 'lambda', 'scores' or 'both'")
    }
    if (is.null(name)) {
        name = 1
    }
    else if ((is.character(name) && length(name) == 1)) {
        if (length(which(names(SFD) == name)) == 1) {
            name = which(names(SFD) == name)
        }
        else if (length(which(names(SFD) == name)) == 0) {
            stop(paste(name, "doesn't not exists. Change name for an existing nameable name."))
        }
        else if (length(which(names(SFD) == name)) == 0) {
            stop("There are more than one nameable with the same name")
        }
    }
    if ((is.null(name) || !(is.numeric(name) && length(name) ==
        1))) {
        stop("Wrong class of name object")
    }
    if (is.null(fill.all)) {
        fill.all = TRUE
    }
    else if (!((isTRUE(fill.all) || isFALSE(fill.all)) && length(fill.all) ==
        1)) {
        stop("Wrong class of fill.all object")
    }
    if (!(inherits(model, "variogramModel") || inherits(model,
        "list"))) {
        stop("Wrong class of model, model should be of class variogramModel or a list of them (use vgm of gstat package) ")
    }
    else if (inherits(model, "list") && !all(sapply(model, inherits,
        "variogramModel"))) {
        stop("Wrong class of model, each element of list should be of class variogramModel (use vgm of gstat package)")
    }
    else if (inherits(model, "list") && (length(model) != ncol(as.data.frame(SFD[[name]]$fpca$scores)))) {
        stop("length of list of models must be equal to number of harmonics of the choosen variable ")
    }
    else if (inherits(model, "variogramModel") && !(fill.all ||
        (ncol(as.data.frame(SFD[[name]]$fpca$scores)) == 1))) {
        stop("If model is not a list and there are more than one nharm of that variable, then fill.all must be TRUE or you can create a list of models with the same number of harmonics")
    }
    puntaje = SFD[[name]]$fpca$scores
    rownames(puntaje) = SFD[[name]]$coordsnames
    puntajes = as.data.frame(puntaje)
    coordinates(puntajes) = SFD[[name]]$coords
    if (fill.all == T) {
        v = list()
        fv = list()
        for (i in 1:ncol(puntajes)) {
            v[[i]] = variogram(puntajes[[i]] ~ 1, puntajes)
            fv[[i]] = fit.variogram(v[[i]], model[[i]], fit.method = 6)
        }
    }
    else {
        v = list()
        fv = list()
        for (i in 1:ncol(puntajes)) {
            v[[i]] = variogram(puntajes[[i]] ~ 1, puntajes)
            fv[[i]] = fit.variogram(v[[i]], model[[i]], fit.method = 6)
        }
    }
    if (method == "lambda" || method == "both") {
        matdis = as.matrix(dist(SFD[[name]]$coords))
        matdis_pred = as.matrix(dist(rbind(SFD[[name]]$coords,
            newcoords)))[(nrow(matdis) + 1):(nrow(matdis) + nrow(newcoords)),
            1:nrow(matdis)]
        if (nrow(newcoords) == 1) {
            matdis_pred <- t(as.matrix(matdis_pred))
        }
        omegas <- list()
        omega = omega <- matrix(0, nrow = nrow(matdis), ncol = ncol(matdis))
        for (i in 1:ncol(puntajes)) {
            omegas[[i]] = variogramLine(model[[i]], dist_vector = matdis,
                covariance = T)
            omega = omega + omegas[[i]]
        }
        vectores_c <- list()
        vector_c = matrix(0, ncol = nrow(matdis), nrow = nrow(matdis_pred))
        for (i in 1:ncol(puntajes)) {
            vectores_c[[i]] = variogramLine(model[[i]], dist_vector = matdis_pred,
                covariance = T)
            vector_c = vector_c + vectores_c[[i]]
        }
        lambda <- solve(omega) \%*\% t(vector_c)
        lambda_data = as.data.frame(lambda)
        colnames(lambda_data) = rownames(newcoords)
        rownames(lambda_data) = SFD[[name]]$coordsnames
        lambda_var <- list()
        for (i in 1:ncol(puntajes)) {
            lambda_var[[i]] <- SFD[[name]]$fpca$values[i] - 2 *
                (vectores_c[[i]] \%*\% lambda) + t(lambda) \%*\%
                (omegas[[i]]) \%*\% lambda
        }
        lambda_var_data <- matrix(0, nrow = nrow(newcoords),
            ncol = ncol(puntajes))
        for (i in 1:ncol(puntajes)) {
            for (j in 1:nrow(newcoords)) {
                lambda_var_data[j, i] <- lambda_var[[i]][j, j]
            }
        }
        lambda_var_data = data.frame(lambda_var_data, VTotal = rowSums(lambda_var_data))
        rownames(lambda_var_data) = rownames(newcoords)
        out_lambda = list(lambda_pred = lambda_data, lambda_varpred = lambda_var_data)
        class(out_lambda) = "lambda_pred"
    }
    if (method == "scores" || method == "both") {
        colnames(newcoords) = c("x", "y")
        coordinates(newcoords) = ~x + y
        K = list()
        for (i in 1:ncol(puntajes)) {
            K[[i]] <- gstat::krige(puntajes[[i]] ~ 1, puntajes,
                newcoords, model = model[[i]], beta = 0)
        }
        pred = K[[1]]$var1.pred
        if (ncol(puntajes) > 1) {
            for (i in 2:ncol(puntajes)) {
                pred = cbind(pred, K[[i]]$var1.pred)
            }
        }
        pred = as.data.frame(pred)
        colnames(pred)[1] = "V1"
        rownames(pred) = rownames(newcoords)
        varpred = K[[1]]$var1.var
        if (ncol(puntajes) > 1) {
            for (i in 2:ncol(puntajes)) {
                varpred = cbind(varpred, K[[i]]$var1.var)
            }
        }
        varpred = data.frame(varpred, VTotal = rowSums(varpred))
        colnames(varpred)[1] = "V1"
        rownames(varpred) = rownames(newcoords)
        out_scores = list(scores_pred = pred, scores_varpred = varpred)
        class(out_scores) = "scores_pred"
    }
    if (method == "both") {
        out <- list(SFD = SFD, KS_scores = out_scores, KS_lambda = out_lambda,
            model = model, name = name)
        class(out) = "KS_pred"
    }
    else if (method == "lambda") {
        out <- list(SFD = SFD, KS_lambda = out_lambda, model = model,
            name = name)
        class(out) = "KS_pred"
    }
    else {
        out <- list(SFD = SFD, KS_scores = out_scores, model = model,
            name = name)
        class(out) = "KS_pred"
    }
    return(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{~Spatial functional prediction}
% \keyword{~Linear combination}
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

