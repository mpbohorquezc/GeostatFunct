ggmap_KS <-
function(KS, map_path=NULL, window_time = NULL, method = "lambda", map_n = 5000, zmin = NULL, zmax = NULL, graph = "plotly"){
  if (!is.null(map_path)){
    if(is.character(map_path)){
      map <- rgdal::readOGR(map_path)
    }else{
      map = map_path
    }
  }else{
    mx <- min(KS$SFD[[1]]$coords[,1])
    Mx <- max(KS$SFD[[1]]$coords[,1])
    my <- min(KS$SFD[[1]]$coords[,2])
    My <- max(KS$SFD[[1]]$coords[,2])
    map <- sp::SpatialPolygons(list(
    sp::Polygons(list(
    sp::Polygon(matrix(c(mx,my,Mx,my,Mx,My,mx,My),byrow = T,ncol = 2))),1)
    ))
  }

  newcoords <- sp::spsample(map, map_n, type = "regular")
  newcoords <- as.data.frame(newcoords)
  colnames(newcoords) <- colnames(KS$SFD[[1]]$coords)

  KS_SFD <- KS_scores_lambdas(KS$SFD, newcoords, model = KS$model, method = method, name = KS$name)

  SFD <- recons_fd(KS_SFD)

  if(is.null(window_time)) {
    times <- SFD$basis$rangeval[1]
  } else if (!(all(window_time >= SFD$basis$rangeval[1]) && all(window_time <= SFD$basis$rangeval[2]))) {
    stop(paste("window_time is out of bounds: Must be some value(s) between ", SFD$basis$rangeval[1], "and ",SFD$basis$rangeval[2]))
  }  else {
    times <- sort(window_time)
  }

  eval <- fda::eval.fd(times, SFD)

  melt_s <- suppressWarnings(reshape::melt(eval))

  melt_s$X2 <- as.factor(melt_s$X2)

  melt_s$X <- as.factor(melt_s$X2)
  levels(melt_s$X) <- newcoords[,1]

  melt_s$Y <- as.factor(melt_s$X2)
  levels(melt_s$Y) <- newcoords[,2]

  names(melt_s) = c("Time","Prediction","Value", "X", "Y")

  melt_s$Time <- as.factor(melt_s$Time)

  graf <- list()
  if(is.null(zmin)){zmin = min(melt_s$Value)}
  if(is.null(zmax)){zmax = max(melt_s$Value)}

  for(i in 1:nlevels(melt_s$Time)){

    melt_s_2 <- melt_s[melt_s$Time == i,]

    if (graph == 'plotly'){
    graf[[i]] <- dplyr::`%>%`(plotly::plot_ly(
      x = as.numeric(as.character(melt_s_2$X)), #melt_s_2$X,
      y = as.numeric(as.character(melt_s_2$Y)), #melt_s_2$Y,
      z = melt_s_2$Value,
      type = "heatmap",
      colorbar = list(title = "Prediction"),
      reversescale = T,
      zmin = zmin,
      zmax = zmax
    ),
      plotly::layout(
        title = paste("Prediction - Time = ", times[i]),
        xaxis = list(showticklabels = FALSE), yaxis = list(showticklabels = FALSE),
        scene = list(aspectration = list(x = 1, y = 1))
      ))
    }
    if (graph == 'gg'){
      graf[[i]] <- ggplot2::ggplot(data = NULL,
                                   aes(x = as.numeric(as.character(melt_s_2$X)), #melt_s_2$X,
                                       y = as.numeric(as.character(melt_s_2$Y))))+ #melt_s_2$Y,
        ggplot2::geom_tile(aes(fill = melt_s_2$Value))+
        ggplot2::labs(fill = "Prediction",title = paste("Prediction - Time = ", times[i]),
                      x = '',y = '',color = NULL,lwd = NULL)+
        ggplot2::scale_fill_viridis_c(direction = -1,limits = c(zmin,zmax)) +
        ggplot2::coord_fixed() +
        ggplot2::theme(plot.background = element_blank(),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       panel.border = element_blank(),
                       axis.line = element_blank(),
                       axis.text = element_blank())
    }

  }

  return(graf)

}
