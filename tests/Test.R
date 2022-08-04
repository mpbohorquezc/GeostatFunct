library(gstat)
library(sp)
library(fda)
library(reshape)
#library(geoR)
library(ggplot2)

# Cargar Datos
PM10 = read.table("PartMat10.txt",head=T,dec=",")
coord = read.table("Coordenadas_planas_Estaciones_Dama_2014.txt",dec=",",sep="\t",header=T)

# Recibir los datos, suavizarlos y ACP
SFD_PM10 <- SpatFD(PM10, coords = coord[, -1], basis = "Bsplines", nbasis = 91, lambda = 0.00002, nharm = 2)

# Modelar los semivariogramas para cada comp. principal
puntaj = data.frame(scores(SFD_PM10)[[1]])

# puntaj1_s = as.geodata(puntaj, coords.col = 1:2, data.col = 3)
# puntaj2_s = as.geodata(puntaj, coords.col = 1:2, data.col = 4)
#
# f1var=variog(puntaj1_s, pairs.min = 4)
# f2var=variog(puntaj2_s, pairs.min = 4)
#
# plot(f1var,main="variog f1")
# plot(f2var,main="variog f2")

#ini1 = eyefit(f1var)
#in12 = eyefit(f2var)

# Definición de los modelos
modelos <- list(vgm(psill = 2634000, "Exp", range = 2103.25, nugget =  0),
                vgm(psill = 101494.96, "Exp", range = 1484.57, nugget = 0))


# Genera los scores y los lambdas para predecir en nuevas coordenadas

#method = "lambda"
KS_SFD_PM10_l <- KS_scores_lambdas(SFD_PM10, coord[1:2,-1]+100 , model = modelos)
#method = "scores"
KS_SFD_PM10_sc <- KS_scores_lambdas(SFD_PM10, coord[1:2,-1]+100 , method = "scores", model = modelos)
#method = "both"
KS_SFD_PM10_both <- KS_scores_lambdas(SFD_PM10, coord[1:2,-1]+100 , method = "both", model = modelos)


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

# Reconstruye los objetos funcionales a partir de los scores y lambdas

#method = "lambda"
curves_PM10_l <- recons_fd(KS_SFD_PM10_l)

#method = "scores"
curves_PM10_sc <- recons_fd(KS_SFD_PM10_sc)
#method = "both"
curves_PM10_both <- recons_fd(KS_SFD_PM10_both)
plot(curves_PM10_both$fd_scores)

# Plot resultados

#method = "lambda"
plot(curves_PM10_l)
ggplot_KS(KS_SFD_PM10_l)
ggplot_KS(KS_SFD_PM10_l, show.varpred = F) #Ocultar la varianza de la predicción

#method = "scores"
plot(curves_PM10_sc)
ggplot_KS(KS_SFD_PM10_sc)

#method = "both"
ggplot_KS(KS_SFD_PM10_both) #Imprime dos plots, uno con la predicción de cada método
ggplot_KS(KS_SFD_PM10_both,
          main = "Plot 1 - Using Scores",
          main2 = "Plot 2 - Using Lambda",
          ylab = "PM10")


# MAPAS
library(rgdal)
library(sf)
library(plotly)
library(sp)

ggmap_KS(KS_SFD_PM10_l,
         map_path = "BOGOTA/Bogota.shp",
         window_time = c(1, 5, 9, 13, 17, 21),
         zmin = 25,
         zmax = 100)


ggmap_KS(KS_SFD_PM10_l,
         map_path = "BOGOTA/Bogota.shp",
         window_time = c(7108),
         method = "lambda",
         zmin = 50,
         zmax = 120)

ggmap_KS(KS_SFD_PM10_l,
         map_path = "BOGOTA/Bogota.shp",
         window_time = c(7108),
         method = "scores",
         zmin = 50,
         zmax = 120)







# COKRIGING

# Cargar Datos Covariable
PM10_co = read.table("PartMat10.txt",head=T,dec=",")

PM10_co = PM10_co + rnorm(10, 7, 2)

coord_co = read.table("Coordenadas_planas_Estaciones_Dama_2014.txt",dec=",",sep="\t",header=T)

# Recibir los datos, suavizarlos, ACP de PM10_co y agregarlo al objeto SFD_PM10
SFD_PM10_co <- SpatFD(PM10_co, coords = coord_co[, -1], basis = "Bsplines", nbasis = 91, lambda = 0.00002, nharm = 2,
                      add = SFD_PM10)


# Extraer los scores de PM10 y de PM10_co
scores(SFD_PM10_co)

# Predicción de PM10_co
KS_SFD_PM10_co  <- KS_scores_lambdas(SFD_PM10_co, coord_co[1:2,-1]+100 , method = "both", model = modelos,
                                     name = "PM10_co")

SFD_PM10_co$name

# Predicción de PM10 usando además la covariable PM10_co #

