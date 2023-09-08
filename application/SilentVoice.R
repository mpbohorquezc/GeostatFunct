# devtools::install_github("mpbohorquezc/GeostatFunct", ref = "main")
library(SpatFD)

library(gstat)
library(sp)
library(rgdal)
library(sf)
library(fda)
library(reshape)
library(tidyverse)

# Load raw data and coordinates
load('SilentVoice_id13.rda')

thinks_a = thinks %>% filter(clase == 'A')
thinks_e = thinks %>% filter(clase == 'E')
thinks_i = thinks %>% filter(clase == 'I')
thinks_o = thinks %>% filter(clase == 'O')
thinks_u = thinks %>% filter(clase == 'U')

# ------------------------------ FUNCTIONAL APPROX ----------------------- # ####
SFD_thinks_a <- SpatFD(thinks_a %>% select(`...1`:`...21`),
                              coords = coord, basis = "Bsplines", nbasis = 91,
                              lambda = 0.00002, nharm = 2,name = 'A')
SFD_thinks_e <- SpatFD(thinks_e %>% select(`...1`:`...21`),
                              coords = coord, basis = "Bsplines", nbasis = 91,
                              lambda = 0.00002, nharm = 2,name = 'E')
SFD_thinks_i <- SpatFD(thinks_i %>% select(`...1`:`...21`),
                              coords = coord, basis = "Bsplines", nbasis = 91,
                              lambda = 0.00002, nharm = 2,name = 'I')
SFD_thinks_o <- SpatFD(thinks_o %>% select(`...1`:`...21`),
                              coords = coord, basis = "Bsplines", nbasis = 91,
                              lambda = 0.00002, nharm = 2,name = 'O')
SFD_thinks_u <- SpatFD(thinks_u %>% select(`...1`:`...21`),
                            coords = coord, basis = "Bsplines", nbasis = 91,
                            lambda = 0.00002, nharm = 2,name = 'U')

ptj_a = data.frame(scores(SFD_thinks_a)[[1]])
ptj_e = data.frame(scores(SFD_thinks_e)[[1]])
ptj_i = data.frame(scores(SFD_thinks_i)[[1]])
ptj_o = data.frame(scores(SFD_thinks_o)[[1]])
ptj_u = data.frame(scores(SFD_thinks_u)[[1]])

sp::coordinates(ptj_a) = sp::coordinates(ptj_e) =
  sp::coordinates(ptj_i) = sp::coordinates(ptj_o) =
  sp::coordinates(ptj_u) = ~x+y

# --------------------------- VARIOGRAMS MODELS -------------------------- # ####

# an example is done, the models are loaded
f1var <- gstat::variogram(sc1_A~1,data = ptj_a,cutoff = 80)
f2var <- gstat::variogram(sc2_A~1,data = ptj_a,cutoff = 48)

f1m = gstat::vgm(psill = 10, "Per", range = 90, nugget =  1,
                 add.to = gstat::vgm(psill = 50,'Gau',range = 20,nugget = 1))
plot(f1var,f1m)

f2m = gstat::vgm(psill = 397, "Gau", range = 148, nugget = 0)
plot(f2var,f2m)

# models_a <- list(f1m,f2m)

load('Models_id13.rda')

# ---------------------- KRIGING AND IMAGES PRODUCTION ------------------- # ####
KS_SFD_thinks_a <- KS_scores_lambdas(SFD_thinks_a, coord ,
                                          method = "scores", model = models_a)
KS_SFD_thinks_e <- KS_scores_lambdas(SFD_thinks_e, coord ,
                                          method = "scores", model = models_e)
KS_SFD_thinks_i <- KS_scores_lambdas(SFD_thinks_i, coord ,
                                          method = "scores", model = models_i)
KS_SFD_thinks_o <- KS_scores_lambdas(SFD_thinks_o, coord ,
                                          method = "scores", model = models_o)
KS_SFD_thinks_u <- KS_scores_lambdas(SFD_thinks_u, coord ,
                                          method = "scores", model = models_u)
N_img = 1000
t = seq(SFD_thinks_a[[1]]$data_fd$basis$rangeval[1],
        SFD_thinks_a[[1]]$data_fd$basis$rangeval[2],length.out = N_img)


theme_set(cowplot::theme_nothing())
tema = theme(legend.position = 'none',
             plot.margin = unit(c(0,0,0,0), "null"),
             panel.background = element_rect(fill = "transparent", colour = NA),
             plot.background = element_blank(),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             panel.border = element_blank(),
             axis.line = element_blank(),
             axis.text = element_blank())

dir.create(dirname(paste0('./A/')), recursive=TRUE)
dir.create(dirname(paste0('./E/')), recursive=TRUE)
dir.create(dirname(paste0('./I/')), recursive=TRUE)
dir.create(dirname(paste0('./O/')), recursive=TRUE)
dir.create(dirname(paste0('./U/')), recursive=TRUE)
for (j in 1:N_img){
  ggsave(plot = ggmap_KS(KS_SFD_thinks_a,
                         map_path = NULL,
                         window_time = t[j],
                         method = "scores",graph = 'gg')[[1]] + ggtitle('') +
           scale_x_continuous(expand=c(0,0)) +
           scale_y_continuous(expand=c(0,0)) +
           tema,filename =paste0('./A/think_a_id13_',j,'.png'),
         width = 400,height = 164,units = 'px',device = 'png',bg = NULL)
  ggsave(plot = ggmap_KS(KS_SFD_thinks_e,
                         map_path = NULL,
                         window_time = t[j],
                         method = "scores",graph = 'gg')[[1]] + ggtitle('') +
           scale_x_continuous(expand=c(0,0)) +
           scale_y_continuous(expand=c(0,0)) +
           tema,filename =paste0('./E/think_e_id13_',j,'.png'),
         width = 400,height = 164,units = 'px')
  ggsave(plot = ggmap_KS(KS_SFD_thinks_i,
                         map_path = NULL,
                         window_time = t[j],
                         method = "scores",graph = 'gg')[[1]] + ggtitle('') +
           scale_x_continuous(expand=c(0,0)) +
           scale_y_continuous(expand=c(0,0)) +
           tema,filename =paste0('./I/think_i_id13_',j,'.png'),
         width = 400,height = 164,units = 'px')
  ggsave(plot = ggmap_KS(KS_SFD_thinks_o,
                         map_path = NULL,
                         window_time = t[j],
                         method = "scores",graph = 'gg')[[1]] + ggtitle('') +
           scale_x_continuous(expand=c(0,0)) +
           scale_y_continuous(expand=c(0,0)) +
           tema,filename =paste0('./O/think_o_id13_',j,'.png'),
         width = 400,height = 164,units = 'px')
  ggsave(plot = ggmap_KS(KS_SFD_thinks_u,
                         map_path = NULL,
                         window_time = t[j],
                         method = "scores",graph = 'gg')[[1]] + ggtitle('') +
           scale_x_continuous(expand=c(0,0)) +
           scale_y_continuous(expand=c(0,0)) +
           tema,filename =paste0('./U/think_u_id13_',j,'.png'),
         width = 400,height = 164,units = 'px')
  print(j)
}
