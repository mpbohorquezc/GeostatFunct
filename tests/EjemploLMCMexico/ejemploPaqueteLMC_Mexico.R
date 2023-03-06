rm(list=ls())

library(sp)
library(geoR)
library(gstat)
# library(sgeostat)
# library(geostatsp)
library(fda)
# library(fda.usc)
# library(splines)
# library(colorspace)
# library(dplyr)

# Datos
PM10=read.table("tests/EjemploLMCMexico/PM10.txt",head=T,dec=".")
RAMA_Coordenadas=read.table("tests/EjemploLMCMexico/RAMA_PM10_coordenadas.txt",head=T,row.names=1,dec=",")
estaciones=colnames(PM10)
CoordenadasPM10=RAMA_Coordenadas[estaciones,]
MPM10=as.matrix(PM10,nrow=nrow(PM10),ncol=18,dimnames=c(rownames(PM10[1,4344]),colnames=colnames(PM10)))
PM10spat=SpatFD(MPM10,CoordenadasPM10,basis="Bsplines",nbasis=21,lambda=0.00002,nharm = 1)
scoresPM10spat=scores(PM10spat)
newcoorden=data.frame(x=c(500000,490000),y=c(2000000,1990000))
matrix(c(500000,2000000),nrow=1,ncol=2,byrow=T)

NO2=read.table("tests/EjemploLMCMexico/NO2.txt",head=T,dec=".")
estacionesNO2=read.table("tests/EjemploLMCMexico/estacionesNO2.txt",head=T,dec=".")
colnames(NO2)=rownames(estacionesNO2)
MNO2=as.matrix(NO2,nrow=4292,ncol=13,dimnames=c(rownames(NO2),colnames=colnames(NO2)))
PM10_MNO2_spat=SpatFD(MNO2, coords = estacionesNO2, basis = "Bsplines", nbasis = 27, lambda = 0.00002, nharm = 4,add=PM10spat)
score_PM10_MNO2=scores(PM10_MNO2_spat)
score_PM10_MNO2[[1]]
score_PM10_MNO2[[2]]

coordinates(newcoorden)
PM10_scores_geodata=as.geodata(scoresPM10spat[[1]],data.col=6)
#First eigenvalue 412808.0. This is me maximum for the sill of the variogram of the first score vector.
#Second eigenvalue 3.139178e+04. This is me maximum for the sill of the variogram of the second score #vector.
v1=variog(PM10_scores_geodata,trend="1st",max.dist = 38000)
#x11()
#plot(v1)
#eyefit(v1)
##model
modelPM10_scores1=list(vgm(412808,"Gau",12698.19))
class(modelPM10_scores1)
##Second dimension
PM10_scores_geodata2=as.geodata(score_PM10_MNO2[[1]],data.col=6)
v2=variog(PM10_scores_geodata2,trend="1st",max.dist = 38000)
#x11()
#plot(v2)
#eyefit(v2)
modelPM10_scores2=vgm(31392, "Ste", 18558.89)
model_ind_PM10_sc1_sc2=list(modelPM10_scores1,modelPM10_scores2)
# ERROR
KS_scores_lambdas(PM10spat, newcoords=newcoorden, model=modelPM10_scores1, method = "scores", fill.all = NULL)

#datosf=fdata(MPM10,argvals=1:nrow(MPM10))
nbasis <-21
hourange <- c(1,nrow(MPM10))
lambda=0.000001
harmaccelLfd <- vec2Lfd(c(1,10), hourange)
hourbasis_Bsplines <- create.bspline.basis(hourange,nbasis)
PM10_fdPar_Bspline<-fdPar(fdobj=hourbasis_Bsplines,Lfdobj=harmaccelLfd,lambda)
PM10_fd_Bspline <- smooth.basis(argvals=1:nrow(MPM10),MPM10,PM10_fdPar_Bspline)
PM10_fd_Bspl=PM10_fd_Bspline$fd
PM10_FPC_bspl=pca.fd(PM10_fd_Bspl,centerfns=T)
#plot.pca.fd(PM10_FPC_bspl)
puntajePM10=PM10_FPC_bspl$scores
colnames(puntajePM10)=c("SC_PM10_1","SC_PM10_2")
rownames(puntajePM10)=estaciones
puntajesPM10=as.data.frame(puntajePM10)

puntajesPM10_1=data.frame(CoordenadasPM10[,1:2],puntajesPM10)
#write.table(puntajesPM10_1,"puntajesPM10.txt",dec = ".",row.names = TRUE,col.names = TRUE)
puntajesg_SC_PM10_1=as.geodata(puntajesPM10_1)
#puntajesg_SC_PM10_2=as.geodata(puntajesPM10_1,data.col=4)
plot(variog(puntajesg_SC_PM10_1))
#plot(variog(puntajesg_SC_PM10_2))
#plot(variog(puntajesg_SC_PM10_1,max.dist=50000))

nbasis <-27
hourange <- c(1,nrow(MNO2))
lambda=0.000001
harmaccelLfd <- vec2Lfd(c(1,10), hourange)
hourbasis_Bsplines <- create.bspline.basis(hourange,nbasis)
NO2_fdPar_Bspline<-fdPar(fdobj=hourbasis_Bsplines,Lfdobj=harmaccelLfd,lambda)
NO2_fd_Bspline <- smooth.basis(argvals=1:nrow(MNO2),MNO2,NO2_fdPar_Bspline)
NO2_fd_Bspl=NO2_fd_Bspline$fd
NO2_FPC_bspl=pca.fd(NO2_fd_Bspl,centerfns=T)
#plot.pca.fd(NO2_FPC_bspl)
puntaje=NO2_FPC_bspl$scores
colnames(puntaje)=c("SC_NO2_1","SC_NO2_2")
rownames(puntaje)=colnames(NO2)
puntajesNO2=as.data.frame(puntaje)
puntajesNO2=data.frame(estacionesNO2[,1:2],puntajesNO2)
puntajesg_SC_NO2_1=as.geodata(puntajesNO2) #79.6%
variPuntaje1_NO2=variog(puntajesg_SC_NO2_1)
plot(variPuntaje1_NO2)
variPuntaje1_NO2=variog(puntajesg_SC_NO2_1,max.dist=60000)
#puntajesg_SC_NO2_2=as.geodata(puntajesNO2,data.col=4) #8.8%
#variPuntaje2_NO2=variog(puntajesg_SC_NO2_2,max.dist=60000)
#plot(variog(puntajesg_SC_NO2_2,max.dist=60000,trend="1st"))
#plot(variPuntaje2_NO2)

coordinates(puntajesPM10_1)=puntajesPM10_1[c("X", "Y")]
coordinates(puntajesNO2)=puntajesNO2[c("X", "Y")]

# ERROR
puntajes=data.frame(RAMA_Coordenadas[,1:2],puntaje)
SC_PM10_1.vgm = variogram(SC_PM10_1~1,puntajes)
SC_PM10_2.vgm = variogram(SC_PM10_2~1,puntajes)
#cross variogram
g = gstat(NULL,"SC_PM10_1",SC_PM10_1~1,puntajesPM10_1,model=vgm(885093,"Exp",23317.05))
g = gstat(g,"SC_NO2_1",SC_NO2_1~1,puntajesNO2,model=vgm(311722.8,"Wav",9408.63))
v = variogram(g)
plot(v,g$model)

################
################
################
################
################

g = gstat(g, fill.all = TRUE)
g.fit = fit.lmc(v, g)
g.fit
plot(v,g.fit)


curve(52476.917-gaussiano(x,52476.917,32258.17),0	,57000,xlab=expression(italic(h)),ylab="",cex.lab=1.5,cex.axis=1.2)
legend("bottomright",legend=c(expression(hat(sigma)^2*" = 52476.917"),expression(hat(phi)*"   = 32258.17")),cex=1.2)
loc <- par("usr")
text(loc[1], loc[4], expression(hat(gamma)[bold(italic(f))[italic(2)]]), xpd = T, adj = c(2.2,4.8),cex=2)
puntos=cbind(xx,yy)
xx2=variPuntaje2_NO2$u
yy2=variPuntaje2_NO2$v
puntos=cbind(xx2,yy2)
points(puntos,pch=19)

# write.table(puntajesNO2,"puntajesNO2.txt")
# write.table(puntajesPM10,"puntajesPM10.txt")
puntajesPM10_NO2=puntajesNO2[-4,]
puntajesPM10_NO2=puntajesNO2_NO2[-13,]
puntajesPM10_NO2=data.frame(puntajesPM10_NO2,puntajesNO2[,3:4])

coordinates(puntajesPM10_NO2)=puntajesPM10_NO2[c("X", "Y")]
SC_PM10_1.vgm = variogram(SC_PM10_1~1,puntajes)
SC_PM10_2.vgm = variogram(SC_PM10_2~1,puntajes)
#cross variogram
g = gstat(NULL,"SC_PM10_1",SC_PM10_1~1,puntajesPM10_NO2)
g = gstat(g,"SC_PM10_2",SC_PM10_2~1,puntajesPM10_NO2)
g = gstat(g,"SC_NO2_1",SC_NO2_1~1,puntajesPM10_NO2)
g = gstat(g,"SC_NO2_2",SC_NO2_2~1,puntajesPM10_NO2)
v = variogram(g)
g = gstat(g, model = vgm(794543, "Exp",22849.54,0), fill.all = TRUE)
g.fit = fit.lmc(v, g)
g.fit
plot(v,g.fit)

###################GRAFICOS SEMIVARIOGRAMA###########################

genCauchy=function(sigma2,h,phi,k1,k2)
{
sigma2-sigma2*(1+(h/phi)^k2)^(-k1/k2)
}
#score1
variPuntaje1=variog(puntajesg_SC_PM10_1,max.dist=60000)
xx=variPuntaje1$u
yy=variPuntaje1$v
puntos=cbind(xx,yy)

exponencial=function(h,cs,a){cs*exp(-h/a)}
gaussiano=function(h,cs,a){cs*exp(-h^2/a^2)}
curve(genCauchy(786932,x,14433.37,1.23,1.31),0,55000,xlab=expression(italic(h)),ylab="",cex.lab=1.5,cex.axis=1.2)
#title(expression(hat(gamma)[bold(italic(f))[italic(1)]]))
loc <- par("usr")
text(loc[1], loc[4], expression(hat(gamma)[bold(italic(f))[italic(2)]]), xpd = T, adj = c(2.2,4.8),cex=2)
#write.table(puntos,"varigramaSC1PM10MEx.txt")
points(puntos,pch=19)
legend("topleft",legend=c(expression(hat(sigma)^2*" = 786932"),expression(hat(phi)*"   = 14433.37"),expression(hat(kappa)[1]*"  = 1.23"),expression(hat(kappa)[2]*"  = 1.31")),cex=1.1)
#title(ylab=expression(hat(gamma)[bold(italic(f))[italic(2)]]),las=3,cex.lab=1.4,mgp=c(1,0,0))
#par(las=1)
#mtext(expression(hat(gamma)[bold(italic(f))[italic(1)]]), side=2, line=2,cex=1.3)

#write.table(puntajes,"puntajesPM10.txt",dec = ".",row.names = TRUE,col.names = TRUE)
puntajesg_SC_PM10_1=as.geodata(puntajes)
puntajesg_SC_PM10_2=as.geodata(puntajes,data.col=4)
#plot(variog(puntajesg_SC_PM10_1,max.dist=60000))
plot(variog(puntajesg_SC_PM10_2,max.dist=60000))#no se observa autocorrelación espacial en el segundo componente principal
##FPC 1 83.1%
##FPC 1 8%

variPuntaje2=variog(puntajesg_SC_PM10_2,max.dist=60000)
xx2=variPuntaje2$u
yy2=variPuntaje2$v
puntos2=cbind(xx2,yy2)
curve(genCauchy(796586.51,x,9184.87,0.08,2),0,55000,xlab=expression(italic(h)),ylab="",cex.lab=1.5,cex.axis=1.2)
#title(expression(hat(gamma)[bold(italic(f))[italic(1)]]))
loc <- par("usr")
text(loc[1], loc[4], expression(hat(gamma)[bold(italic(f))[italic(2)]]), xpd = T, adj = c(2.2,4.8),cex=2)
#write.table(puntos,"varigramaSC1PM10MEx.txt")
points(puntos2,pch=19)
legend("topleft",legend=c(expression(hat(sigma)^2*" = 786932"),expression(hat(phi)*"   = 14433.37"),expression(hat(kappa)[1]*"  = 1.23"),expression(hat(kappa)[2]*"  = 1.31")),cex=1.1)

curve(106771.07-gaussiano(x,106771.07,25881.81),0,57000,xlab=expression(italic(h)),ylab="",cex.lab=1.5,cex.axis=1.2)
legend("bottomright",legend=c(expression(hat(sigma)^2*" = 106771.07"),expression(hat(phi)*"   = 25881.81")),cex=1.2)
loc <- par("usr")
text(loc[1], loc[4], expression(hat(gamma)[bold(italic(f))[italic(2)]]), xpd = T, adj = c(2.2,4.8),cex=2)

curve(794543.91-exponencial(x,794543.91,22849.54),0,57000,xlab=expression(italic(h)),ylab="",cex.lab=1.5,cex.axis=1.2)
legend("bottomright",legend=c(expression(hat(sigma)^2*" = 794543.91"),expression(hat(phi)*"   = 22849.54")),cex=1.2)
loc <- par("usr")
text(loc[1], loc[4], expression(hat(gamma)[bold(italic(f))[italic(1)]]), xpd = T, adj = c(2.2,4.8),cex=2)

plot(NO2_fd_Bspl)
####Gráfico último semana, parece q aumenta en época de verano
plot(NO2_fd_Bspl,xlim=c(3420,3588),col=rainbow(18),lwd=2,lty=1,ylab="",xlab="hour",cex.lab=1.2,cex.axis=1.2)
#plot(NO2_fd_Bspl,ylab="",xlim=c(4140,4294.5),xlab="hour",cex.lab=1.2,cex.axis=1.2)
lines(NO2_fd_Bspl,col=rainbow(18),lwd=2,lty=1)
title(ylab="NO2  (ppm)",cex.lab=1.2,mgp=c(2.5,0,0))
lines(NO2_fd_Bspl,col=rainbow(9),lwd=2,lty=1)
legend(5995,160,legend=ID_estaciones,lwd=2,cex=1,col=rainbow(9))
#plotfit.fd(MPM10, argvals=1:nrow(MPM10),PM10_fd_Bspl,lwd=1,ylab=" ")
#Para graficar uno solo, por ejemplo la estación 1
#plot(PM10_fd_Bspline$fd[1])
