library(readr)

# Cargar Datos
PM10 <- read.csv("test/EjemploLMCMexico/PM10.txt", sep="")
NO2 <- read.csv("test/EjemploLMCMexico/NO2.txt", sep="")
RAMA_Coordenadas = read.table("test/EjemploLMCMexico/RAMA_PM10_coordenadas.txt", head=T, row.names=1, dec=",")

View(PM10)
View(NO2)
View(RAMA_Coordenadas)

estacionesPM10 = colnames(PM10)
CoordenadasPM10 = RAMA_Coordenadas[estacionesPM10,]
MPM10 = as.matrix(PM10, nrow=nrow(PM10), ncol=ncol(PM10), dimnames=c(rownames(PM10), colnames=colnames(PM10)))

CoordenadasNO2 = RAMA_Coordenadas[c(-4, -14),]
colnames(NO2) = rownames(CoordenadasNO2)
MNO2 = as.matrix(NO2,nrow=nrow(NO2),ncol=ncol(NO2),dimnames=c(rownames(NO2),colnames=colnames(NO2)))

# Recibir los datos, suavizarlos y ACP
SFD_PM10_NO2 <- SpatFD(MPM10, coords = CoordenadasPM10[,1:2], basis = "Fourier", nbasis = 21, lambda = 0.000001, nharm = 2)
SFD_PM10_NO2 <- SpatFD(MNO2, coords = CoordenadasNO2[, 1:2], basis = "Fourier", nbasis = 27, lambda = 0.000001, nharm = 2,
                      add = SFD_PM10_NO2)

# Semivariogramas para cada comp. principal
puntaj_PM10 = data.frame(scores(SFD_PM10_NO2)[[1]])
puntaj_NO2 = data.frame(scores(SFD_PM10_NO2)[[2]])

puntaj_PM10_1 = as.geodata(puntaj_PM10, coords.col = 1:2, data.col = 3)
puntaj_PM10_2 = as.geodata(puntaj_PM10, coords.col = 1:2, data.col = 4)

puntaj_NO2_1 = as.geodata(puntaj_NO2, coords.col = 1:2, data.col = 3)
puntaj_NO2_2 = as.geodata(puntaj_NO2, coords.col = 1:2, data.col = 4)


PM10_f1var=variog(puntaj_PM10_1)
PM10_f2var=variog(puntaj_PM10_2)

NO2_f1var=variog(puntaj_NO2_1)
NO2_f2var=variog(puntaj_NO2_2)

par(mfrow = c(2,2))

plot(PM10_f1var,main="variog PM10_f1")
plot(PM10_f2var,main="variog PM10_f2")
plot(NO2_f1var,main="variog NO2_f1")
plot(NO2_f2var,main="variog NO2_f2")

# 2 Componentes por cada variable
coordinates(puntaj_PM10) <- ~ X + Y
coordinates(puntaj_NO2) <- ~ X + Y

# Valor propio comp principal
#PM10
SFD_PM10_NO2$MPM10$fpca$varprop
SFD_PM10_NO2$MPM10$fpca$values
#NO2
SFD_PM10_NO2$MNO2$fpca$varprop
SFD_PM10_NO2$MNO2$fpca$values

# Semivariograma
gst = gstat(NULL,"PM10_1", sc1_MPM10~1, data = puntaj_PM10)
gst = gstat(gst,"PM10_2", sc2_MPM10~1, data = puntaj_PM10 )
gst = gstat(gst,"NO2_1", sc1_MNO2~1, data = puntaj_NO2)
gst = gstat(gst,"NO2_2", sc2_MNO2~1, data = puntaj_NO2)

v = variogram(gst)
plot(v)

# Modelo Semivariograma gstat
model1 <- gstat::vgm(647677.1,"Gau",23317.05)
model1 <- gstat::vgm(127633,"Wav",9408.63, add.to = model1)
model1

gfit <- fit.lmc(v, gst, model = model1, correct.diagonal=1.01,
                fit.method = 6,
                fit.ranges = F)
gfit
plot(v, gfit$model)

# Predict
newcoords <- matrix(c(509926, 2179149), nrow = 1)
colnames(newcoords) = c('x','y')
newcoords = as.data.frame(newcoords)
coordinates(newcoords) = ~x+y
z = predict(gfit, newdata = newcoords)
z


# Matriz de distancias

distan11 = as.matrix(dist(CoordenadasPM10))
distan22 = as.matrix(dist(CoordenadasNO2))

CoordenadasPM10_NO2 = rbind(CoordenadasPM10, CoordenadasNO2)
nrow(CoordenadasPM10)
nrow(CoordenadasNO2)

distan12 = as.matrix(dist(CoordenadasPM10_NO2))
distan12 = distan12[(nrow(CoordenadasPM10)+1):nrow(CoordenadasPM10_NO2),
                    1:nrow(CoordenadasPM10)]

# Matriz de corregionalización

# "Gau" range = 23317.05

B1 <- matrix(c(gfit$model$PM10_1[1,2], gfit$model$PM10_1.PM10_2[1,2], gfit$model$PM10_1.NO2_1[1,2], gfit$model$PM10_1.NO2_2[1,2],
  0, gfit$model$PM10_2[1,2], gfit$model$PM10_2.NO2_1[1,2], gfit$model$PM10_2.NO2_2[1,2],
  0, 0, gfit$model$NO2_1[1,2], gfit$model$NO2_1.NO2_2[1,2],
  0, 0, 0,gfit$model$NO2_2[1,2]), nrow = 4)
library(gdata)
upperTriangle(B1) <- lowerTriangle(B1)
B1

# "Wav" range = 9408.63
B2 <- matrix(c(gfit$model$PM10_1[2,2], gfit$model$PM10_1.PM10_2[2,2], gfit$model$PM10_1.NO2_1[2,2], gfit$model$PM10_1.NO2_2[2,2],
               0, gfit$model$PM10_2[2,2], gfit$model$PM10_2.NO2_1[2,2], gfit$model$PM10_2.NO2_2[2,2],
               0, 0, gfit$model$NO2_1[2,2], gfit$model$NO2_1.NO2_2[2,2],
               0, 0, 0,gfit$model$NO2_2[2,2]), nrow = 4)
upperTriangle(B2) <- lowerTriangle(B2)
B2








# xxxxxxx

modeloscok <- list(vgm(psill = 381766.8, "Gau", range = 7125.97, nugget =  0),
                vgm(psill = 41295.88, "Gau", range = 7125.97, nugget = 0),
                vgm(psill = 213668.1, "Gau", range = 7113.069, nugget = 0),
                vgm(psill = 78562.56, "Gau", range = 7113.069, nugget = 0))

# Matriz de corregionalización

B1=matrix(c(0.61,0.3,0.4,0.3,
            0.3,0.001,0.1,0,
            0.4,0.1,0.05,0.1,
            0.3,0,0.1,0.006),byrow=T,nrow=4)
B1=nearPD(B1)$mat

B2=matrix(c(0.003,-0.1,-0.16,-0.62,
            -0.1,0.24,0.12,0.4,
            -0.16,0.12,0,0,
            -0.62,0.4,0,0.002),byrow=T,nrow=4)
B2=nearPD(B2)$mat

B3=matrix(c(0.001,0,-0.17,0.35,
            0,0.002,-0.2,0,
            -0.17,-0.2,0.19,0.10,
            0.35,0,0.10,0.001),byrow=T,nrow=4)

B3=nearPD(B3)$mat
B4=matrix(c(0.05,0,0.22,0.7,
            0,0,0.07,0,
            0.22,0.07,0,0.23,
            0.7,0,0.23,0.57),byrow=T,nrow=4)
B4=nearPD(B4)$mat

B4=matrix(c(0.0,0,0.0,0.0,
            0,0,0.0,0,
            0,0,0.57,0.23,
            0,0,0.23,0.57),byrow=T,nrow=4)
B4=nearPD(B4)$mat

pa = matrix(c(1,2,3,1,2,3), nrow = 3, ncol = 3)
pb = matrix(c(1,2,2,1,2,2), nrow = 3, ncol = 3)
pa*pb

# Argumentos #
#SpatFD   SFD_PM10_NO2
#Modelos  modeloscok
#Matrices B1 B2 B3 B4



PM10_1=function(B1,B2,B3,B4,x){B1[1,1]*Modelo_qc_1(x)+B2[1,1]*Modelo_qc_2(x)+B3[1,1]*Modelo_u2_1(x)+B4[1,1]*Modelo_fs_1(x)}
PM10_2=function(B1,B2,B3,B4,x){B1[2,2]*Modelo_qc_1(x)+B2[2,2]*Modelo_qc_2(x)+B3[2,2]*Modelo_u2_1(x)+B4[2,2]*Modelo_fs_1(x)}
NO2_1=function(B1,B2,B3,B4,x){B1[3,3]*Modelo_qc_1(x)+B2[3,3]*Modelo_qc_2(x)+B3[3,3]*Modelo_u2_1(x)+B4[3,3]*Modelo_fs_1(x)}
NO2_2=function(B1,B2,B3,B4,x){B1[4,4]*Modelo_qc_1(x)+B2[4,4]*Modelo_qc_2(x)+B3[4,4]*Modelo_u2_1(x)+B4[4,4]*Modelo_fs_1(x)}
PM10_1_PM10_2=function(B1,B2,B3,B4,x){B1[1,2]*Modelo_qc_1(x)+B2[1,2]*Modelo_qc_2(x)+B3[1,2]*Modelo_u2_1(x)+B4[1,2]*Modelo_fs_1(x)}
PM10_1_NO2_1=function(B1,B2,B3,B4,x){B1[1,3]*Modelo_qc_1(x)+B2[1,3]*Modelo_qc_2(x)+B3[1,3]*Modelo_u2_1(x)+B4[1,3]*Modelo_fs_1(x)}
PM10_1_NO2_2=function(B1,B2,B3,B4,x){B1[1,4]*Modelo_qc_1(x)+B2[1,4]*Modelo_qc_2(x)+B3[1,4]*Modelo_u2_1(x)+B4[1,4]*Modelo_fs_1(x)}
PM10_2_NO2_1=function(B1,B2,B3,B4,x){B1[2,3]*Modelo_qc_1(x)+B2[2,3]*Modelo_qc_2(x)+B3[2,3]*Modelo_u2_1(x)+B4[2,3]*Modelo_fs_1(x)}
PM10_2_NO2_2=function(B1,B2,B3,B4,x){B1[2,4]*Modelo_qc_1(x)+B2[2,4]*Modelo_qc_2(x)+B3[2,4]*Modelo_u2_1(x)+B4[2,4]*Modelo_fs_1(x)}
NO2_1_NO2_2=function(B1,B2,B3,B4,x){B1[3,4]*Modelo_qc_1(x)+B2[3,4]*Modelo_qc_2(x)+B3[3,4]*Modelo_u2_1(x)+B4[3,4]*Modelo_fs_1(x)}







g = gstat(NULL,"PM10_1", X1~1, data = puntaj_PM10, model=vgm(885093,"Exp",23317.05))
g = gstat(g,"NO2_1", X1~1, data = puntaj_NO2, model=vgm(311722.8,"Wav",9408.63))
#g = gstat(g, c("PM10_1","NO2_1"), data = puntaj_NO2, model=vgm(180722.8,"Exp",5408.63))

modelos2 <- vgm(885093,"Exp",23317.05, add.to = vgm(311722.8,"Wav",9408.63))
gfit <- fit.lmc(v, g, modelos2, fit.method = 6, correct.diagonal=1.01)
v = variogram(gfit)
plot(v,
     model = gfit$model,
     pl = T,
     xlab = "Distancias",
     ylab = "Semivarianza")

newcoords <- matrix(c(509926, 2179149), nrow = 1)
colnames(newcoords) = c('x','y')
newcoords = as.data.frame(newcoords)
coordinates(newcoords) = ~x+y
z = predict(gfit, newdata = newcoords)




gfit <- fit.lmc(v, g,
                fit.method = 2,
                fit.ranges = F)

gfit
matdis = as.matrix(dist(RAMA_Coordenadas[1:2]))
variogramLine(gfit, matdis)



sum(SFD_PM10_NO2$MPM10$fpca$harmonics$coefs[,1] %*% t(SFD_PM10_NO2$MNO2$fpca$harmonics$coefs[,1]))
sum(SFD_PM10_NO2$MPM10$fpca$harmonics$coefs[,2] %*% t(SFD_PM10_NO2$MNO2$fpca$harmonics$coefs[,1]))
sum(SFD_PM10_NO2$MPM10$fpca$harmonics$coefs[,2] %*% t(SFD_PM10_NO2$MNO2$fpca$harmonics$coefs[,2]))

sum(SFD_PM10_NO2$MPM10$fpca$harmonics$coefs[,1] %*% t(SFD_PM10_NO2$MPM10$fpca$harmonics$coefs[,1]))
sum(SFD_PM10_NO2$MPM10$fpca$harmonics$coefs[,2] %*% t(SFD_PM10_NO2$MPM10$fpca$harmonics$coefs[,2]))

SFD_PM10_NO2$MPM10$data_fd$basis

SFD_PM10_NO2$MPM10$fpca$harmonics * SFD_PM10_NO2$MPM10$fpca$harmonics
plot(SFD_PM10_NO2$MPM10$fpca$harmonics, plot = FALSE)

SFD_PM10_NO2$MPM10$fpca$harmonics

str(SFD_PM10_NO2$MNO2$fpca$harmonics)
SFD_PM10_NO2$MPM10$data_fd$coefs

sum(eval.fd(seq(0,1,by = 0.1), SFD_PM10_NO2$MPM10$fpca$harmonics[1]) * eval.fd(seq(0,1,by = 0.1), SFD_PM10_NO2$MPM10$fpca$harmonics[1]))

sum((eval.basis(seq(1,4344,by = 0.1), SFD_PM10_NO2$MPM10$fpca$harmonics[1]$basis)) * (eval.basis(seq(1,4344,by = 0.1), SFD_PM10_NO2$MPM10$fpca$harmonics[1]$basis)))

eval1 = eval.basis(seq(1,4344,by = 0.1), SFD_PM10_NO2$MPM10$fpca$harmonics[1]$basis)
eval2 = eval.basis(seq(1,4344,by = 0.1), SFD_PM10_NO2$MPM10$fpca$harmonics[2]$basis)
sum(matplot(eval2, type = "l")*matplot(eval2, type = "l"))
sum(eval1*eval2)
