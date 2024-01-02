#Construye una grilla para x
#ini_x: punto de x donde inicia la grilla
#fin_x: punto de x donde finaliza la grilla
#nx: número de puntos discretizados en que se divide la grilla para x
#nsp: Número de puntos espaciales del proceso
#spmod: model espacial que explica la variablidad del proceso
#mu: media del proceso
#sigma2: Varianza
#phi: 

library(mvtnorm)
library(geoR)
library(Matrix)

grilla1 <- function(iniz=0, finz=1, nz=100)
#Genera la grilla para una variable
{
    z=seq(iniz, finz, length.out=nz)
    return(z)
}


grilla <- function(ReadXY=TRUE, XY, ini_x=0, fin_x=1, nx=10, ini_y=0, fin_y=1, ny=10)
{
 if (ReadXY)
 {
    x<-XY[,1]
    y<-XY[,2]
    nx <- length(x)
    ny <- length(y)
 }
 else
 {
   x <- grilla1(ini_x, fin_x, nx)
   y <- grilla1(ini_y, fin_y, ny)
 }
   nsp <- nx * ny
      
   grilla_1=expand.grid(x,y)
   return(grilla_1)
}


XYcoord <- function(ReadXY=TRUE, XY=matrix(0,nrow=3, ncol=3), ini_x=0, fin_x=1, nx=10, ini_y=0, fin_y=1, ny=10)
#This function generates the distances matrix from a grid of spatial coordinates. The grid can be provide using a matrix or simulated.
#ReadXY: Logical value to indicate if the spatial coordinates are readed from a matrix or they are simulated
#X:       n x 2 Matrix with the coordinates of n sites. It is obligatory if ReadXY is TRUE. The j_th row represents the spatial
#          coordinates (x,y) for the point S_j 

{
   ##############################################################
   #                      Spatial coordinates                   #
   ##############################################################   
   if (ReadXY)
   {
      cat("Se usará la información de las localizaciones espaciales suministrada en la matrix X","\n")
      if (!is.matrix(XY)) stop(" !!! The COORDINATE MATRIX must be a MATRIX !!! ")
      if (length(which (is.nan(XY))) > 0  )     stop(" !!! The COORDINATES MATRIX can not contain NAN values !!! ")
      if (ncol(XY) !=2) stop(" The COORDINATE MATRIX  must have 2 columns")      
      grilla_1 <- XY
   }
   else
   {
      cat("Las localizaciones espaciales se generarán discretizando los intervalos para x e y","\n")      
      if(fin_x <= ini_x || fin_y <= ini_y ) 
               stop(" !!! poorly defined grid parameters. The end point must be greater than the start point !!! ")
      if(nx <= 0 || ny <= 0)
            stop(" !!! the number of points for the grid must be positive !!!")
      grilla_1 <- grilla(ReadXY, XY, ini_x, fin_x, nx, ini_y, fin_y, ny)
   }
   
   
   nsp <- nrow(grilla_1)
   distancias_1=as.matrix(dist(grilla_1))
   return(list(nsp=nsp, grilla_1=grilla_1, distMat=distancias_1))
}


SpCovMat <- function(ReadCovMat=FALSE, B=matrix(0, nrow=3, ncol=3), distMat, Spmod="exponential", Sigma2=1, Phi=5)
#This function READ the spatial covariance matrix if ReadCovMat is TRUE, else it is estimated.
#distMat:    Distance Matrix between the sites                       
#ReadCovMat: Boolean variable to determine if the spatial covariance matrix is read from a matriz given (ReadCovMat=TRUE)
#            it is estimated from a list of spatial models (ReadCovMat=FALSE)). 
#            If the ReadCovMat=TRUE, spatial covariance matrix (B) must be provided. 
#            Else, the spatial covariance matriz is estimated and the followinf information of the P
#            spatial processes must be provided.            
#                        Spmod:   Name of the spatial model associated 
#                        Sigma2:  Sill of the spatial model associated 
#                        Phi:     Range of the spatial model associated 

{                                               
   if (ReadCovMat)
   {
      ##########################################################################
      #                     TO READ THE COVARIANCE MATRIX                      #
      ##########################################################################   
       if (!is.matrix(B))
          stop ("¡¡¡ ERROR. The provided matrix IS NOT COVARIANCE MATRIX !!!")
       if (nrow(B) != ncol(B))
          stop ("¡¡¡ ERROR. THE MATRIX IS NOT SQUARED. !!!") 
       if  (length(which(diag(B) < 0)) > 0)
          stop ("¡¡¡ ERROR. THERE ARE VARIANCE values NON POSITIVE !!!")                   
       sigmaX <- B
   }
   else
   {

       sigmaX <- cov.spatial(distMat,cov.model=Spmod,cov.pars=c(Sigma2, Phi))
    }
   #En caso de ser no definida positiva se aproxima por la matrix def + más cercana
   sigmaX <- as.matrix(nearPD(sigmaX)$mat)
   return(sigmaX)
}


SimSpatScoreP <- function(ReadXY=TRUE,  XY, ini_x=0, fin_x=1, nx=10, ini_y=0, fin_y=1, ny=10, 
                          ReadCovMat=FALSE, CovLST, Spmod=rep("exponential", 3), Mu=rep(0, 3), Sigma2=rep(1,3), Phi=rep(5,3))
#ARGUMENTOS ASOCIADOS A LA GRILLA ESPACIAL.                          
#           ReadXY:  Logical value to indicate if the data are readed from a matrix or generated.
#                    TRUE: The spatial grid is provided by X matrix
#                    FALSE: The spatial grid is generated using ini_x, fin_x, nx, ini_y, fin_y y ny.
#                    The Following arguments IS obligatory if ReadXY=TRUE: 
#                        XY:       n x 2 Matrix with the coordinates of n sites. It is obligatory if ReadXY is TRUE.
#                    The following arguments are obligatory if ReadXY=FALSE: 
#                        ini_x:   Inititial value for the grid in the x axis.
#                        fin_x:   End value for the grid in the x axis.
#                        nx:      Number of points in the x axis.
#                        ini_y:   Inititial value for the grid in the y axis.
#                        fin_y:   End value for the grid in the y axis.
#                        ny:      Number of points in the y axis.
#ARGUMENTOS ASOCIADOS A LA FORMA DE SUMINISTRAR LA MATRIZ DE COVARIANZA
#           ReadCovMat: Logical Value to indicate if the covariance matrices are readed from a LIST of matrices (one for eah score)
#                       or are generated using spatial models, one for each score.        
#                     TRUE: Indicates that the spatial dependence structure for each score is provided using a LIST of matrices.
#                           The matrices are associated with the scores in the order that are reported
#                     FALSE: The covariance matrix is generated from spatial models. 
#                            The models and their parameters are asociated with the scores according to
#                            the order in that the models are reported.
#ARGUMENTOS ASOCIADOS A LOS MODELOS QUE GENERAN LA MATRIZ DE COVARIANZA
#           Spmod:   Vector with the names of the P spatial models associated a each of the P scores.
#           Mu:      Vector with the means of the Scores associated to P scores
#           Sigma2:  Vector with the variances of the P spatial models associated to each score.  
#           Phi:     Vector with the range of the P spatial models associated to each score.  
#NOTAS:
#If the user needs to use the spatial grid stored in X, ReadXY must be TRUE
#If the user needs to discretize the spatial information for the x and y intervals, ReadXY is FALSE.

{
   Grid     <- XYcoord(ReadXY, XY, ini_x, fin_x, nx, ini_y, fin_y, ny)
   nsp      <- Grid$nsp
   grilla_1 <- Grid$grilla_1
   distMat  <- Grid$distMat
   
   if (ReadCovMat==TRUE)   
   {
      P <- length(CovLST)
      for (p in 1:P)
          if (nrow(CovLST[[p]]) != nsp) 
             stop("!!! ERROR.THE NUMBER OF ROWS AND COLUMNS OF THE COVARIANCE MATRICES MUST BE EQUAL TO THE NUMBER OF SPATIAL POINTS")
   }
   else
   {
      ##############################################################################
      #                     THE COVARIANCE MATRICES ARE ESTIMATED                  #
      ##############################################################################   
      if ((length(Mu) != length(Spmod)) || (length(Mu) != length(Sigma2)) || (length(Mu) != length(Phi)))
         stop("!!! ERROR. THE DIMENSIONS OF THE PARAMETERS ASSOCIATED WITH THE SPATIAL PROCESS ARE DIFFERENT.¡¡¡")
      if ( length(which((Sigma2 <= 0))) > 0 )
         stop("!!! ERROR. SOME VARIANCE values ARE NON POSITIVE ¡¡¡")
      if ( length(which((Phi <= 0))) > 0 )
         stop("!!! ERROR. SOME RANGE VALUES ARE NON POSITIVE .¡¡¡")         
      P <- length(Mu)
   }   
   
   Z        <- matrix(0, nrow=nsp, ncol=P)

   for (p in 1:P)
   {            
      if (ReadCovMat==TRUE)
          CovMat <- SpCovMat(ReadCovMat, CovLST[[p]])
      else
          CovMat <- SpCovMat(distMat=distMat, Spmod=Spmod[p], Sigma2=Sigma2[p], Phi=Phi[p])
          
      Z[,p] <- rmvnorm(1, mean=rep(Mu[p],nsp), sigma=CovMat)
   }
   sim1  <- cbind(grilla_1,Z)
   name1 <- rep(" ",P)
   for (p in 1:P)
       name1[p] <- paste("Score",p)
   names(sim1)=c("x","y",name1)    
   return(sim1)
}


#EXAMPLES 

#EXAMPLE 1.
#ReadXY=FALSE: The Spatial Coordinates are taken from a grid and the distance matrix is calculated from this grid. 
#    By default, for the x and y axis, 10 points are generated between 0 and 1, therefore ther are 100 sites.
#    These values can be changed using ini_, fin_x, nx and ini_y, fin_y and ny.
#ReadCovMat=FALSE: The Covariance Matrices for the scores are generated using spatial models. 
#    By default, the number of scores is 3. The spatial model assoctiated for each score is exponential with mean=0, sill=1 y range=5.
#**********************************************************************************************************************
sC <- SimSpatScoreP(ReadXY=FALSE, ReadCovMat=FALSE)


#EXAMPLE 2.
#ReadXY=TRUE: The Spatial Coordinates are taken from a MATRIX (n x 2), where n is the number of sites.
#ReadCovMat=FALSE: The Covariance Matrices for the scores are generated using spatial models. 
#    By default these models are exponential with mean=0,  sill=1 y range=5.
#**********************************************************************************************************************
c1 <- c(0,0,0,1,1,1)
c2 <- c(0,1,2,0,1,2)
XY <- cbind(c1,c2)
sC <- SimSpatScoreP(ReadXY=TRUE, XY=XY, ReadCovMat=FALSE)


#EXAMPLE 3.
#ReadXY=TRUE: The Spatial Coordinates are taken from a MATRIX (n x 2), where n is the number of sites.
#ReadCovMat=TRUE: The Covariance Matrices for the scores are provided USING A LIST OF MATRICES
#**********************************************************************************************************************
c1 <- c(0,1)
c2 <- c(1,2)
XY <- cbind(c1,c2)
CovLst <- list(matrix(1,nrow=4,ncol=4), matrix(2,nrow=4,ncol=4))
sC <- SimSpatScoreP(ReadXY=TRUE, XY=XY, ReadCovMat=TRUE, CovLST=CovLst)


#EXAMPLE 4.
#ReadXY=FALSE: The Spatial Coordinates are taken from a grid and the distance matrix is calculated from this grid. 
#    By default, for the x and y axis, 10 points are generated between 0 and 1, therefore ther are 100 sites.
#    These values can be changed using ini_, fin_x, nx and ini_y, fin_y and ny.
#ReadCovMat=TRUE: The Covariance Matrices for the scores are provided USING A LIST OF MATRICES
#**********************************************************************************************************************
CovLst <- list(matrix(1,nrow=4,ncol=4), matrix(2,nrow=4,ncol=4))
sC <- SimSpatScoreP(ReadXY=FALSE, nx=2, ny=2, ReadCovMat=TRUE, CovLST=CovLst)

#EXAMPLE 5.
c1 <- c(0,1)
c2 <- c(1,2)
XY <- cbind(c1,c2)
#sC <- SimSpatScoreP(ReadXY=FALSE, ReadCovMat=FALSE, Mu<-c()
