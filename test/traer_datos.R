library(readxl)
library(dplyr)
#library(proj4)


traer_datos <- function(variables,
                        fecha,
                        hora) {

    estaciones <- read.csv("Data_frame/estaciones_planas.csv", header = T)
    rutas <- paste("Data_frame/20RAMA/2020",
                    variables,
                    ".xls",
                    sep = "")

    xlsx <- lapply(rutas, read_excel)
    xlsx <- lapply(xlsx, as.data.frame)
    xlsx <- lapply(xlsx, mutate, FECHA = as.character(FECHA))
    xlsx <- lapply(xlsx, filter, FECHA == fecha, HORA == hora)
    xlsx <- sapply(xlsx, select, -FECHA, -HORA, simplify = T)
    frame <- as.data.frame(t(bind_rows(xlsx)))
    colnames(frame) <- variables
    frame["cve_estac"] <- row.names(frame)

    frame_final <- estaciones %>%
                    select(cve_estac, X, Y, alt) %>%
                    inner_join(frame) %>%
                    rename(Estacion = cve_estac, altitud = alt)

    frame_final[frame_final == -99.0] <- NA

    frame_final <- frame_final[!rowSums(frame_final[variables],
                    na.rm = T) == 0.0, ]

    return(frame_final)

    }

# read.csv("Data_frame/estaciones_planas.csv", he)
# Acá tengo que preguntar lo de las coordenadas porque hay algunas
# estaciones que no estná las poryecciones.
# project(cbind(datos$longitud, datos$latitud),
#         proj = "+proj=utm +zone=34 +south +ellps=WGS84
#         +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
# 36      SAC -99.00938 19.34561    2293 0.1  1   6   NA    9   0 45   6   NA