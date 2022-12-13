# test
#setwd("/home/samuel/Documentos/U/GstatFD/GeostatFunct-main/tests")


source("../R/SpatFD.R")
source("../R/OSD_scores_lambda.R")

library(rgdal)
library(gstat)
library(sp)

# bogota_shp <- rgdal::readOGR("BOGOTA/Bogota.shp")

# save(bogota_shp,file = "Bogota.rda")

data(Bogota)

vgm_model  <- gstat::vgm(psill = 5.665312,
                  model = "Exc",
                  range = 8000,
                  kappa = 1.62,
                  add.to = vgm(psill = 0.893,
                               model = "Nug",
                               range = 0,
                               kappa = 0))

vgm_model_2 <- list(gstat::vgm(5.66,"Sph",8000),
                    gstat::vgm(5.66,"Exp",8000))

my.CRS <- sp::CRS("+init=epsg:21899") # https://epsg.io/21899

bogota_shp <- sp::spTransform(bogota_shp,my.CRS)
target <- sp::spsample(bogota_shp,n = 100, type = "random") # The set of points in which we want to predict optimally.
old_stations <- sp::spsample(bogota_shp,n = 3, type = "random") # The set of stations that are already fixed.
my_grid <- sp::spsample(bogota_shp,n = 3000, type = "regular")


# Case 1 ####

FD_optimal_design(k = 10, s0 = target,model = vgm_model,#method = "scores",
               map = bogota_shp,plt = T) -> res1
res1
class(res1$new_stations)
class(res1$fixed_stations)
plot(res1)

# Case 2 ####

FD_optimal_design(k = 10, s0 = target,model = vgm_model,
               map = bogota_shp,plt = T,#method = "scores",
               fixed_stations = old_stations) -> res2
res2
class(res2$fixed_stations)
plot(res2)

# Case 3 ####

FD_optimal_design(k = 10, s0 = target,model = vgm_model,
               grid = my_grid,plt = T,nharm = 3,
               fixed_stations = old_stations) -> res3

res3
plot(res3)

# Case 4 ####

FD_optimal_design(k = 4, s0 = target, model = vgm_model_2,
                  map = bogota_shp,method = "scores", plt = T) -> res4
res4
plot(res4)

