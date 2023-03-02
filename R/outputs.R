\name{outputs}
\alias{outputs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Summary of SpatFD objects
%%  ~~function to do ... ~~
}
\description{
This functions shows a summary of the main objects of SpatFD objects.
}
\usage{
summary(X)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     A SpatFD object
}
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
\value{For each variable included in the SpatFd object, this functions return: Head of data, Coordinates, Eigenvalues, Mean coefficients, Proportion of explained variance by each component
%%
%%
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Joan Nicolás Castro
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
PM10=read.table("PM10.txt",head=T,dec=".")
RAMA_Coordenadas=read.table("RAMA_PM10_coordenadas.txt",head=T,row.names=1,dec=",")
estaciones=colnames(PM10)
CoordenadasPM10=RAMA_Coordenadas[estaciones,]
MPM10=as.matrix(PM10,nrow=nrow(PM10),ncol=18,dimnames=c(rownames(PM10[1,4344]),colnames=colnames(PM10)))
PM10spat=SpatFD(MPM10,CoordenadasPM10,basis="Bsplines",nbasis=21,lambda=0.00002,nharm=3)
summary(PM10spat)

summary <- function(SpatFD){
  if (class(SpatFD)=="SpatFD"){
    i=1
    for (i in 1:length(SpatFD)){
      var = SpatFD[[i]]$variable_name
      df = SpatFD[[i]]$data
      coor = SpatFD[[i]]$coords
      ev = SpatFD[[i]]$fpca$values
      meanfd = SpatFD[[i]]$fpca$meanfd
      varprop = SpatFD[[i]]$fpca$varprop

      cat("# ",var,"\n")

      cat("## Data","\n")
      print(rbind(head(df, 4),
                  rep("...", times=dim(df)[2]),
                  tail(df, 4)))

      cat("\n","## Coordinates","\n")
      print(rbind(head(coor,4),
                  rep("...", times=dim(coor)[2])))

      cat("\n","## Eigenvalues","\n")
      print(rbind(head(data.frame(ev),4),
                  "..."))

      cat("\n","## Mean coefficients","\n")
      print(rbind(head(data.frame(meanfd$coefs), 4),
                  "...",
                  tail(data.frame(meanfd$coefs),4)))

    cat("\n","## Proportion of explained variance by component","\n")
    print(rbind(head(data.frame(varprop))))
    cat("\n","\n")

    i=i+1
  }}
  if (class(SpatFD)=="KS_pred"){
    if (is.null(SpatFD$KS_scores)&&!is.null(SpatFD$KS_lambda)){
      lambda_pred = SpatFD$KS_lambda$lambda_pred
      lambda_varpred = SpatFD$KS_lambda$lambda_varpred
      model = SpatFD$model

      cat("\n","## Lambda values","\n")
      print(lambda_pred)

      cat("\n","## Lambda var_predicted","\n")
      print(lambda_varpred)

      cat("\n","## Models","\n")
      for (i in 1:length(model)){
        cat("The model ",i,"is: \n")
        print(model[[i]])}
    }
    if (!is.null(SpatFD$KS_scores)&&is.null(SpatFD$KS_lambda)){
      scores_pred = SpatFD$KS_scores$scores_pred
      scores_varpred = SpatFD$KS_scores$scores_varpred
      model = SpatFD$model

      cat("\n","## Scores","\n")
      print(scores_pred)

      cat("\n","## Scores var_predicted","\n")
      print(scores_varpred)

      cat("\n","## Models","\n")
      for (i in 1:length(model)){
        cat("The model ",i,"is: \n")
        print(model[[i]])}

    }
    if (!is.null(SpatFD$KS_scores)&&!is.null(SpatFD$KS_lambda)){
      scores_pred = SpatFD$KS_scores$scores_pred
      scores_varpred = SpatFD$KS_scores$scores_varpred
      lambda_pred = SpatFD$KS_lambda$lambda_pred
      lambda_varpred = SpatFD$KS_lambda$lambda_varpred
      model = SpatFD$model

      cat("\n","## Lambda values","\n")
      print(lambda_pred)

      cat("\n","## Lambda var_predicted","\n")
      print(lambda_varpred)

      cat("\n","## Scores","\n")
      print(scores_pred)

      cat("\n","## Scores var_predicted","\n")
      print(scores_varpred)

      cat("\n","## Models","\n")
      for (i in 1:length(model)){
        cat("The model ",i,"is: \n")
        print(model[[i]])}
    }}
}

objects <- function(SpatFD){

    data = SFD_PM10$PM10$data
    coords = SFD_PM10$PM10$coords
    data_fd = SFD_PM10$PM10$data_fd
    coefs = SFD_PM10$PM10$data_fd$coefs
    basis = SFD_PM10$PM10$data_fd$basis

    cat("## Data","\n")
    print(rbind(head(data, 3),
                rep("...", times=dim(data)[2]),
                tail(data, 3)))

    cat("\n","## Coords","\n")
    print(rbind(head(coords, 3),
                rep("...", times=dim(coords)[2]),
                tail(coords, 3)))

    cat("\n","## Data fd","\n")
    print(data_fd)

    cat("\n","## Coefs","\n")
    print(rbind(head(coefs, 3),
                rep("...", times=dim(coords)[2]),
                tail(coefs, 3)))

    cat("\n","## Basis","\n")
    print(head(basis))
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
