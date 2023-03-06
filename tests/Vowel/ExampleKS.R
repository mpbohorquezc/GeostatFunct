library(gstat)
library(geoR)

#Read data and coordinates
coordins=read.table("coordin.txt",head=T,dec=",")
subject1=read.table("subjet1Epsd.txt",head=T,dec=",")
subject1_dat1_a=subject1[1:228,]
save(coordins, subject1,  file = "Subject1_vowels.RData")
save(coordins, subject1,  file = "Subject1_vowels.rda")
#load("Subject1_vowels.rda")
# Build SpatFD object. This function build curves and perform FPCA
SFD_sub1 <- SpatFD(subject1_dat1_a, coords = coordins, basis = "Bsplines", nbasis = 15, lambda = 0.00002, nharm = 2)
Scores_sub1 = data.frame(scores(SFD_sub1)[[1]])

##########################################################################################
#Create variogram object using gstat package, if models and parameters are already known.
##### If don't at the end there are some guide for this task.
##############################################################
vgm_scores1=vgm(107.57,"Exc", 15.43,kappa = 1.35)
vgm_scores2=vgm(43.56,"Exp",29.09)
var_models <- list(vgm_scores1,vgm_scores2)

# Genera los scores y los lambdas para predecir en nuevas coordenadas
x=seq(16,64,by=10.5)
y=seq(0,48,by=10.5)
newcoords=expand.grid(x,y)
colnames(newcoords)=c("x","y")
#coordinates(newcoords)=c("x","y")

#method = "lambda"
KS_SFD_sub1_l <- KS_scores_lambdas(SFD_sub1, newcoords, model = var_models)
#method = "scores"
KS_SFD_sub1_sc <- KS_scores_lambdas(SFD_sub1, newcoords, method = "scores", model = var_models)
#method = "both"
KS_SFD_sub1_both <- KS_scores_lambdas(SFD_sub1, newcoords , method = "both", model = var_models)

# Retorna a SFD_PM10, scores, lambdas y las varianzas de predicción

#method = "lambda"
KS_SFD_sub1_l$KS_lambda$lambda_pred
KS_SFD_sub1_l$KS_lambda$lambda_varpred

#method = "scores"
KS_SFD_sub1_sc$KS_scores$scores_pred
KS_SFD_sub1_sc$KS_scores$scores_varpred

#method = "both"

KS_SFD_sub1_both$KS_lambda$lambda_pred
KS_SFD_sub1_both$KS_lambda$lambda_varpred
KS_SFD_sub1_both$KS_scores$scores_pred
KS_SFD_sub1_both$KS_scores$scores_varpred

# Reconstruye los objetos funcionales a partir de los scores y lambdas

#method = "lambda"
curves_sub1_l <- recons_fd(KS_SFD_sub1_l)
plot(curves_sub1_l)
#method = "scores"
curves_sub1_sc <- recons_fd(KS_SFD_sub1_sc)
plot(curves_sub1_sc)
#method = "both"
curves_sub1_both <- recons_fd(KS_SFD_sub1_both)

ggplot_KS(KS_SFD_sub1_both) #Imprime dos plots, uno con la predicción de cada método

# Plot resultados

#method = "lambda"

ggplot_KS(KS_SFD_sub1_l)
ggplot_KS(KS_SFD_sub1_sc)

ggplot_KS(KS_SFD_sub1_l, show.varpred = F) #Ocultar la varianza de la predicción

#method = "scores"

ggplot_KS(KS_SFD_sub1_sc)

###method = "both"

ggplot_KS(KS_SFD_sub1_sc) #Imprime dos plots, uno con la predicción de cada método,
ggplot_KS(KS_SFD_sub1_both,
          main = "Plot 1 - Using Scores",
          main2 = "Plot 2 - Using Lambda",
          ylab = "Subject1")

###########
ggmap_KS(KS_SFD_sub1_l,map_path = NULL,
         window_time = c(1, 5, 9, 13, 17, 21),
         zmin = min(subject1_dat1_a),
         zmax = max(subject1_dat1_a))

ggmap_KS(KS_SFD_sub1_l,map_path = NULL,
         window_time = c(10),
         method = "lambda",
         zmin = min(subject1_dat1_a),
         zmax = max(subject1_dat1_a))

ggmap_KS(KS_SFD_sub1_l,map_path = NULL,
         window_time = c(150),
         method = "scores",
         zmin = min(subject1_dat1_a),
         zmax = max(subject1_dat1_a))

####################################?????????????????????????????????
#If you still do not know these models, you can use geoR or gstat packages to find them.
#Finding semivariogram models for spatial score vector of each dimension included.

######################################################################################
#Fitting Semivariogram model for 1st score, 1st subject
######################################################################################

Scores_sub1 = data.frame(scores(SFD_sub1)[[1]])
Scores_sub1_score1_geo=as.geodata(Scores_sub1,data.col = 3)
plot(Scores_sub1_score1_geo)
plot(Scores_sub1_score1_geo,trend="1st")
vari_score1=variog(Scores_sub1_score1_geo,trend="1st",max.dist=66)
x11()
plot(vari_score1)
ini=eyefit(vari_score1)
#1:228
#cov.model sigmasq   phi tausq kappa kappa2   practicalRange
#1 powered.exponential  107.57 15.43     0  1.35   <NA> 34.7801345375644
#In this case, the fitted model in geoR was powered exponential (stable), and now we convert it to a vgm object, which is the object-class for variograms in gstat package.

######################################################################################
#Fitting Semivariogram model for 2nd score
######################################################################################
Scores_sub1_score2_geo=as.geodata(Scores_sub1,data.col = 4)
plot(Scores_sub1_score2_geo)
model_sub1_score2=lm(X2~x+y,Scores_sub1)
Scores_sub1$resi_score2=residuals(modelscore2)
Scores_sub1_score2_geo=as.geodata(Scores_sub1,data.col = 5)
plot(Scores_sub1_score2_geo)
vari_score2=variog(Scores_sub1_score2_geo)
x11()
plot(vari_score2)
ini=eyefit(vari_score2)
ini
# ini
#cov.model sigmasq   phi tausq kappa kappa2   practicalRange
#1 exponential   43.56 29.09     0  <NA>   <NA> 87.1458518364096

######################################################################################################
###Defining first parameters for build functions: type of basis, number of basis, smoothing parameter
######################################################################################################

################################################################################
###From discrete observations to functional data.  Build and plot using Fourier
################################################################################

Msubject1_dat1_a=as.matrix(subject1_dat1_a,nrow=228,ncol=21,dimnames=c(rownames(subject1_dat1_a),colnames=colnames(subject1_dat1_a)))
nbasis <-13
hourange <- c(1,nrow(Msubject1_dat1_a))
lambda=0.00002
harmaccelLfd <- vec2Lfd(c(1,21), hourange)
hourbasis <- create.fourier.basis(hourange,period=228/1,nbasis)
sub1_a_fdPar_Fou <- fdPar(fdobj=hourbasis,Lfdobj=harmaccelLfd,lambda)
sub1_a_fd <- smooth.basis(Msubject1_dat1_a,argvals=1:nrow(Msubject1_dat1_a),sub1_a_fdPar_Fou)
sub1_a_fd_Fou= sub1_a_fd$fd
plot(sub1_a_fd_Fou)
lines(sub1_a_fd_Fou,col=heat.colors(21),lwd=2,lty=1)
legend("topright",colnames(subject1),lwd=2,col=heat.colors(21),cex=0.6)
x11()15
plotfit.fd(Msubject1_dat1_a, argvals=1:nrow(Msubject1_dat1_a),sub1_a_fd_Fou,lwd=1,ylab=" ")

###############################################################################
#From discrete observations to functional data.  Build and plot using Bsplines
###############################################################################

Msubject1_dat1_a=as.matrix(subject1_dat1_a,nrow=228,ncol=21,dimnames=c(rownames(subject1_dat1_a),colnames=colnames(subject1_dat1_a)))
nbasis <-19
hourange <- c(1,nrow(Msubject1_dat1_a))
lambda=0.000001
harmaccelLfd <- vec2Lfd(c(1,21), hourange)
hourbasis_Bsplines <- create.bspline.basis(hourange,nbasis)
subject1_dat1_a_fdPar_Bspline<-fdPar(fdobj=hourbasis_Bsplines,Lfdobj=harmaccelLfd,lambda)
subject1_dat1_a_fd_Bspline <- smooth.basis(argvals=1:nrow(Msubject1_dat1_a),Msubject1_dat1_a,subject1_dat1_a_fdPar_Bspline)
subject1_dat1_a_fd_Bspl=subject1_dat1_a_fd_Bspline$fd
x11()
plot(subject1_dat1_a_fd_Bspl)
lines(subject1_dat1_a_fd_Bspl,col=topo.colors(21),lwd=2,lty=1)
legend("topright",colnames(subject1),lwd=2,col=topo.colors(21),cex=0.6)
plotfit.fd(Msubject1_dat1_a, argvals=1:nrow(Msubject1_dat1_a),subject1_dat1_a_fd_Bspl,lwd=1,ylab=" ")



########################################################################################################################################################
vgm_scores1=vgm(107.57,"Exc", 15.43,kappa = 1.35)
Scores_sub1_g = Scores_sub1
coordinates(Scores_sub1_g)= c("x","y")
Scores_sub1_gstat_score1=gstat(id="X1", formula = X1~x+y,data=Scores_sub1_g,model = vgm_scores1)
v_empir_score1=variogram(Scores_sub1_gstat_score1,cutoff=66)
vari_fitted_score1=fit.variogram(v_empir_score1, vgm_scores1)
attr(vari_fitted_score1, "SSErr")
plot(v_empir_score1,vari_fitted_score1)
#fit.variogram(v_empir_score1, vgm_scores1,fit.sills = T,fit.ranges = T,fit.kappa=T)

vgm_scores2=vgm(43.56,"Exp",29.09)
Scores_sub1_gstat_score2=gstat(id="X2", formula = X2~x+y,data=Scores_sub1_g,model = vgm_scores2)
v_empir_score2=variogram(Scores_sub1_gstat_score2,cutoff=90,boundaries=vari_score2$u)
plot(v_empir_score2)
vari_fitted_score2=fit.variogram(v_empir_score2, vgm_scores2)
attr(vari_fitted_score2, "SSErr")
plot(v_empir_score2,vari_fitted)



###############################################################################
#Exploring variogrqams with gstat
###############################################################################

Scores_sub1_gstat_score1=gstat(id="X1", formula = X1~1,data=Scores_sub1_g,model = vgm(var(Scores_sub1$X1), "Exp", max(dist(Scores_sub1[,1:2]))))

Scores_sub1_gstat_score2=gstat(id="X2", formula = X2~1,data=Scores_sub1_g,model = vgm(var(Scores_sub1$X2), "Exp", max(dist(Scores_sub1[,1:2]))))
#Scores_sub1=gstat(id="X1", formula = X1~1, data=Scores_sub1_g, model = vgm(14687.15, "Exp", 94))
v_empir_score1=variogram(Scores_sub1_gstat_score1,cutoff=max(dist(Scores_sub1[,1:2])))
#View(variogram(Scores_sub1_gstat,cutoff=max(dist(Scores_sub1[,1:2]))))
plot(v_empir_score1,cutoff=max(dist(Scores_sub1[,1:2])))
v_empir_score1=variogram(Scores_sub1_gstat_score1,cutoff=85)
plot(v_empir_score1,cutoff=80)
plot(variogramLine(variomodel_scores1,t="l",ylim=c(0,60000))
     plot(variogramLine(vgm(95000, "Mat",45,kappa=1.4),100),t="l",ylim=c(0,60000))
     points(v_empir_score1$dist,v_empir_score1$gamma)
     plot(variogramLine(vgm(85000, "Gau",105),80),t="l",,ylim=c(0,60000))
     points(v_empir_score1$dist,v_empir_score1$gamma)
     v.fit = fit.variogram(v_empir_score1,variomodel_scores1)
     x11()

