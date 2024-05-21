library(gstat)
library(fda)
library(sp)

# setwd("~/extra/U/GstatFD/Simulación condicional/")


data("CanadianWeather")
dim(CanadianWeather$dailyAv)


# coordenadas
canada.CRS <- CRS("+init=epsg:4608")
coords <- SpatialPoints(CanadianWeather$coordinates,
                        proj4string = CRS("+init=epsg:4326"))
coords <- spTransform(coords,canada.CRS)


# Observaciones
obs <- CanadianWeather$dailyAv[,,1] # Temperatura
dim(obs)

## suavizar
Lfd_obj <- int2Lfd(m = 2)
create.bspline.basis(rangeval = c(1,365),
                     nbasis = 40, norder = 4) -> mi.base

GCV.bsp <- NULL
log_lambdas <- -40:4
for (k in log_lambdas){
  lambda <- exp(k)
  mi.fdPar <- fdPar(mi.base, Lfd_obj, lambda = lambda)
  mi.fd <- smooth.basis(argvals = 1:365,
                        y = obs, fdParobj = mi.fdPar)
  GCV.bsp <- c(GCV.bsp,sum(mi.fd$gcv))
}
log_lam <- log_lambdas[which.min(GCV.bsp)]

mi.fdPar <- fdPar(mi.base, Lfd_obj, lambda = exp(log_lam))
mi.fd <- smooth.basis(argvals = 1:365,
                      y = obs, fdParobj = mi.fdPar)

nbasis <- 5
canada <- mi.fd$fd
canada.pca <- pca.fd(canada,nharm = 10)
base_ort <- canada.pca$harmonics[1:nbasis]
canada_mean <- canada.pca$meanfd
class(canada.pca$scores)

# media

formula2fd <- function(rango, expresion) {
  # Generar puntos equidistantes en el rango dado
  n <- 500  # Número de puntos para evaluar
  x <- seq(rango[1], rango[2], length.out = n)
  
  # Evaluar la expresión en los puntos x_vals
  y_vals <- eval(parse(text = expresion))
  
  # Crear un objeto fd con los valores obtenidos
  basis <- create.bspline.basis(rangeval = rango, nbasis = 30)
  fd_obj <- Data2fd(x, y_vals,basisobj = basis)
  
  return(fd_obj)
}

media <- formula2fd(c(-1,1),"3*sin(x*4)")
# plot(media)

# Simular

source("simulations.R")


## No condicional
vario <- vgm(.25, "Exp", .5, .05)
# plot(vario,cutoff = 2)
nbasis <- 6
sims <- sim_functional_process(10,vario,nbasis,coords,basis = 'Legendre',mu = media)

class(sims)
length(sims)

class(sims[[1]])
# plot(sims[[3]][[1]]$data_fd)


sims <- sim_functional_process(10,vario,nbasis,coords,basis = 'Legendre')

class(sims)
length(sims)

class(sims[[1]])
# plot(sims[[3]][[1]]$data_fd)


# Condicional
vario <- vgm(100, "Exp", 900, 10)
# plot(vario,cutoff = 5000)

new_coords <- spsample(coords,100,type = "regular")
gridded(new_coords) <- TRUE

# par(bty = "l",pty = "s")
# plot(coordinates(coords),pch = '*', col = "blue",
#      xlab = "x",ylab = "y")
# points(coordinates(new_coords),pch = 20, col = "green")
# grid()

length(new_coords)
a <- sim_functional_process(10,vario,nbasis,new_coords,canada,coords)

class(a)
length(a)

x11()
class(a[[1]])
#plot(a[[1]][[1]]$data_fd)


vario <- vgm(100, "Wav", 900, 10)
#plot(vario,cutoff = 5000)
a <- sim_functional_process(10,vario,nbasis,new_coords,canada,coords)

class(a)
length(a)

#x11()
class(a[[1]])
#plot(a[[1]][[1]]$data_fd)


