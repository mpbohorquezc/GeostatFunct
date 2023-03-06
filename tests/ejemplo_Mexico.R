library(gstat)

PM10 = read.table("tests/EjemploLMCMexico/PartMat10.txt",head=T,dec=",")
coord = read.table("tests/EjemploLMCMexico/Coords_DAMABog.txt",dec=",",sep="\t",header=T,row.names=1)

# Create an univariate object using 2 nharm
SFD_PM10 <- SpatFD(PM10, coords = coord, basis = "Bsplines", nbasis = 91, lambda = 0.00002, nharm = 2)
View(SFD_PM10)
summary(SFD_PM10)
objects(SFD_PM10)

# Create a bi-variate object using 3 nharm for the first variable and 1 nharm for the second variable
PM10=read.table("tests/EjemploLMCMexico/PM10.txt",head=T,dec=".")
RAMA_Coordenadas=read.table("tests/EjemploLMCMexico/RAMA_PM10_coordenadas.txt",head=T,row.names=1,dec=",")
estaciones=colnames(PM10)
CoordenadasPM10=RAMA_Coordenadas[estaciones,]
MPM10=as.matrix(PM10,nrow=nrow(PM10),ncol=18,dimnames=c(rownames(PM10[1,4344]),colnames=colnames(PM10)))
PM10spat=SpatFD(MPM10,CoordenadasPM10,basis="Bsplines",nbasis=21,lambda=0.00002,nharm=3)
summary(PM10spat)
objects(PM10spat)

NO2=read.table("tests/EjemploLMCMexico/NO2.txt",head=T,dec=".")
estacionesNO2=read.table("tests/EjemploLMCMexico/estacionesNO2.txt",head=T,dec=".")
nrow(estacionesNO2)
colnames(NO2)=rownames(estacionesNO2)
MNO2=as.matrix(NO2,nrow=4292,ncol=13,dimnames=c(rownames(NO2),colnames=colnames(NO2)))
PM10_MNO2_spatFD=SpatFD(MNO2, coords = estacionesNO2, basis = "Bsplines", nbasis = 27, lambda = 0.00002, nharm = 1,add=PM10spat)
summary(PM10_MNO2_spatFD)
objects(PM10_MNO2_spatFD)
