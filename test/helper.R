source("traer_datos.R")
source("plot_tools.R")

library(sp)
library(gstat)
library(sf)
library(rgdal)
library(ggplot2)
library(plotly)
library(stringr)
library(Matrix)


remove_na <- function(frame, vari_) {

    datos1 <- frame

    bool <- !is.na(datos1@data[vari_])
    datos1@data <- datos1@data[bool, ]
    datos1@coords <- datos1@coords[bool, ]

    return(datos1)

}

square <- function(dimen) {

    return(paste(c(as.character(c(diag(dimen)))), collapse = ","))
}


get_matrix <- function(string, dim) {

    return(matrix(as.numeric(unlist(str_split(string, ","))), nrow = dim))
}

get_matrix("1,0,0,1", 2)

build_vmg <- function(lista_rangos, list_sillas, lista_modelos) {

    vgm <- vgm(list_sillas[1], lista_modelos[1], lista_rangos[1])
    n_mol <- length(lista_modelos)
    if (n_mol > 1) {
        for (ii_model in 2:n_mol) {
            vgm <- vgm(list_sillas[ii_model],
                    lista_modelos[ii_model],
                    lista_rangos[ii_model],
                    add.to = vgm)
            }
    }
    return(vgm)

}

extrac_index <- function(list_pair){

    return(unique(unlist(str_split(list_pair,"-"))))

}



vector_order <- function(total){

    vec <- c(1)
    last <- 1
    for (number in  (total):2) {
        vec <- c(vec, last + number)
        last <- last + number
    }

    jum <- total * (total + 1) / 2
    unique(c(vec, 1:jum))

}

# vector_order(4)

g_model <- function(listas_variables,
                    modelos,
                    matrices,
                    formu,
                    rangos,
                    datos1) {


    datos <- datos1
    coordinates(datos) <- ~ X + Y
    n_variables <- length(listas_variables)

    matrices <- lapply(matrices, get_matrix, dim =  n_variables)
    matrices <- lapply(matrices, function(x){ return(nearPD(x)$mat)})


    cruzados <- matrix(apply(expand.grid(listas_variables, listas_variables), 1,
                        paste, collapse = "-"), nrow = n_variables)



    indices <- lower.tri(matrices[[1]], diag = TRUE)
    cruzados <- cruzados[indices]

    list_id <- sapply(cruzados, extrac_index)

    n_cruzados <- length(cruzados)

    cocurrre_matrix <- sapply(matrices, function(x) {return(x[indices])})

    list_model <- list()
    for (row in 1:n_cruzados) {

        appen <- build_vmg(lista_rangos = rangos,
                list_sillas = cocurrre_matrix[row, ],
                lista_modelos = modelos)

        list_model[[cruzados[row]]] <- appen
    }

    g <- NULL

    a <- 1

    for (cur in vector_order(n_variables)) {

        var <- cruzados[cur]
        if (a < (n_variables + 1)) {

            g <- gstat(g, id = list_id[[var]],
                    formula =  as.formula(formu[a]),
                    data = remove_na(datos, listas_variables[[a]]),
                    model = list_model[[var]])
            a <- a + 1
        } else {
            g <- gstat(g, id = list_id[[var]],
                        model = list_model[[var]])

        }
}


return(g)

}




# g <- g_model(c("O3", "NO2"),
#             c("Nug", "Gau", "Hol"),
#             matrices <- c("1,1,1,1", "25,15,15,25", "25,15,15,25"),
#             formu <- c("O3 ~ 1", "NO2 ~ X + altitud"),
#             rangos <- c(0, 2000, 2000),
#             datos)



# predict(g, newdata = loci)

# plot(variogram(g), model = g$model)

# g <- fit.lmc(variogram(g), g,
#             fit.method = 2,
#             fit.ranges = F)

# plot(variogram(g), model=g$model)



# vgm()
# list_id[[2]]



