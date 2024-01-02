source("traer_datos.R")
source("plot_tools.R")
source("helper.R")

library(sp)
library(gstat)
library(sf)
library(rgdal)
library(ggplot2)
library(plotly)
library(stringr)
library(Matrix)

# Fecha :
# 2020-01-16
# Hora:
# 17
# NO2 ~ X + Y
# O3 ~  Y
# NOX ~  Y

# Sph
# 30,30,10,30,50,10,50,50,35
# Rango 6096

# Hol
# 1,30,20,30,40,40,20,40,20
# 2294

variables <- c("CO", "NO", "NO2",
               "PM10", "PM25",
               "SO2", "O3",
               "NOX", "PMCO")

fecha <- "2020-01-03"
hora <- 20

datos <- traer_datos(variables, fecha, hora)

datos
g <- g_model(c("O3", "NO2", "PM10"),
            c("Nug", "Sph", "Exp"),
            matrices <- c("1,1,1,1,1,1,11,1,1",
                          "1,1,1,1,1,1,11,1,1",
                          "1,1,1,1,1,1,11,1,1"),
            formu <- c("O3 ~ 1", "NO2 ~ X + altitud", "PM10 ~ 1"),
            rangos <- c(0, 10000, 2000),
            datos)

dis_matrix <- as.matrix(dist(cbind(datos$X, datos$Y)))

# Matriz de distancias
# definir un vgm compuesto
vgm_ <- build_vmg(c(0, 10000, 2000),
                  c(1, 1, 1),
                  c("Nug", "Sph", "Exp"))
vgm_

cov_matrix1 <- variogramLine(object = g$model$O3,
                            dist_vector = dis_matrix,
                            covariance = T)
cov_matrix2 <- variogramLine(object = g$model$NO2,
                            dist_vector = dis_matrix,
                            covariance = T)
cov_matrix3 <- variogramLine(object = g$model$PM10,
                            dist_vector = dis_matrix,
                            covariance = T)
cov_matrix4 <- variogramLine(object = g$model$O3.NO2,
                            dist_vector = dis_matrix,
                            covariance = T)
cov_matrix5 <- variogramLine(object = g$model$O3.PM10,
                            dist_vector = dis_matrix,
                            covariance = T)
cov_matrix6 <- variogramLine(object = g$model$NO2.PM10,
                            dist_vector = dis_matrix,
                            covariance = T)

list_sub_matrix <- list(cov_matrix1,
                        cov_matrix2,
                        cov_matrix3,
                        cov_matrix4,
                        cov_matrix5,
                        cov_matrix6)

cov_grade <- matrix(0, ncol = 26 * 5, nrow =  26 * 5)


col_1 <- rbind(cov_matrix1, cov_matrix4, cov_matrix5)
col_2 <- rbind(cov_matrix4, cov_matrix2, cov_matrix6)
col_3 <- rbind(cov_matrix5, cov_matrix6, cov_matrix3)

cov_grande <- cbind(col_1, col_2, col_3)




dim(cov_matrix)
class(g$model)

plot(variogram(g), model = g$model)

g <- fit.lmc(variogram(g), g,
            fit.method = 2,
            fit.ranges = F)


variogram(g), model=g$modelplot()
