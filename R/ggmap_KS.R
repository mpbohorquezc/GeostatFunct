ggmap_KS <-
function(KS, map_path, window_time = NULL, method = "lambda", map_n = 5000, zmin = NULL, zmax = NULL){

  map <- readOGR(map_path)
  newcoords <- sp::spsample(map, map_n, type = "regular")
  newcoords <- as.data.frame(newcoords)
  colnames(newcoords) <- c("X", "Y")

  KS_SFD <- KS_scores_lambdas(KS$SFD, newcoords, model = KS$model, method = method, name = KS$name)

  SFD <- recons_fd(KS_SFD)

  if(is.null(window_time)) {
    times <- SFD$basis$rangeval[1]
  } else if (!(all(window_time >= SFD$basis$rangeval[1]) && all(window_time <= SFD$basis$rangeval[2]))) {
    stop(paste("window_time is out of bounds: Must be some value(s) between ", SFD$basis$rangeval[1], "and ",SFD$basis$rangeval[2]))
  }  else {
    times <- sort(window_time)
  }

  eval <- eval.fd(times, SFD)

  melt_s <- suppressWarnings(melt(eval))

  melt_s$X2 <- as.factor(melt_s$X2)

  melt_s$X <- as.factor(melt_s$X2)
  levels(melt_s$X) <- newcoords$X

  melt_s$Y <- as.factor(melt_s$X2)
  levels(melt_s$Y) <- newcoords$Y

  names(melt_s) = c("Time","Prediction","Value", "X", "Y")

  melt_s$Time <- as.factor(melt_s$Time)

  graf <- list()
  if(is.null(zmin)){zmin = min(melt_s$Value)}
  if(is.null(zmax)){zmax = max(melt_s$Value)}

  for(i in 1:nlevels(melt_s$Time)){

    melt_s_2 <- melt_s[melt_s$Time == i,]

    graf[[i]] <- plot_ly(
      x = as.numeric(as.character(melt_s_2$X)), #melt_s_2$X,
      y = as.numeric(as.character(melt_s_2$Y)), #melt_s_2$Y,
      z = melt_s_2$Value,
      type = "heatmap",
      colorbar = list(title = "Prediction"),
      reversescale = T,
      zmin = zmin,
      zmax = zmax
    ) %>%
      layout(
        title = paste("Prediction - Time = ", times[i]),
        xaxis = list(showticklabels = FALSE), yaxis = list(showticklabels = FALSE),
        scene = list(aspectration = list(x = 1, y = 1))
      )

  }

  return(graf)

}
