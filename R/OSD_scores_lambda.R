.vgm_model.fn <- function(points, s0, model, method_i = "lambda", grid, fixed_stations){
  
  # validation ---------------------------------
  
  # This is an internal function, thus its arguments will be validated in the
  # main function FD_optimal_design.
  
  # Basic objects
  method <- method_i
  N <- nrow(s0)
  points <- matrix(points,ncol=2)
  points <- rbind(fixed_stations,points)
  colnames(points) <- c("x","y")
  # k <- length(points)/2
  k <- nrow(points)
  z <- rep.int(1,k)
  L <- length(model)
  df <- sp::SpatialPointsDataFrame(coords = points,
                                   data = as.data.frame(z))
  targetPoints <- sp::SpatialPoints(s0)
  
  # Replace points by their nearest in the given grid
  D_toGrid <- proxy::dist(points,as.matrix(grid),diag = T)
  closests <- apply(as.matrix(D_toGrid),1,which.min)
  points <- as.matrix(grid)[closests,]
  
  # distance matrix
  dist_matrix <- as.matrix(proxy::dist(points,diag = T))
  dist_s0 <- as.matrix(proxy::dist(points, s0,diag = T))
  
  # ------------------ method lambda
  
  if (method == "lambda")){
    
    # Omega
    omegas <- list()
    omega <- matrix(0, k, k)
    for(j in 1:L){
      omegas[[j]] <- gstat::variogramLine(model[[j]], dist_vector = dist_matrix,coveriance = T)
      omega = omega + omegas[[j]]
      
      # C
      c_vecs <- NULL
      for(l in 1:N)
        c_vecs <- cbind(c_vecs, gstat::variogramLine(model[[j]], dist_vector = dist_s0[,l],coveriance = T))
      
    }
    
    accuracy <- 0
    for(l in 1:N){
      #c_vec <- rep(0,k)
      for(j in 1:L){
        #c_vec <- c_vec + gstat::variogramLine(model[[j]], dist_vector = dist_s0[,l],coveriance = T)
        c_vec <- gstat::variogramLine(model[[j]], dist_vector = dist_s0[,l],coveriance = T)
        accuracy <- c_vec%*%omegas[[j]]%*%t(c_vec)
        # Sugerencia: Ponderar por proporcion de varianza explicada.
      }
    }
    
    return(1/accuracy)
  }
  
  # ------------------ method scores
  
  if (method == "scores")){
    
    # fit variogram
    tot_variance <- 0
    for(j in 1:L){
      
      gstat_obj <- gstat::gstat(
        id = "z",
        formula = z~1,
        model = model[[j]],
        beta = 0,
        data = df)
      
      predictions <- gstat::predict(
        object = gstat_obj,
        newdata = targetPoints
      )
      
      tot_variance <- tot_variance + sum(as.data.frame(predictions)$z.var)
    }
    
    return(tot_variance)
    }
  
}

FD_optimal_design <- function(k,s0,model,fixed_stations = NULL,
                              method = "lambda", grid = NULL,
                              map = NULL,nharm = NULL,plt = F){
  # validation ----------------------------------------------
  
  if (missing(k))
    stop("Missing k.")
  
  if(missing(s0))
    stop("Missing s0.")
  
  if(missing(s0))
    stop("Missing Model.")
  
  if(is.null(grid) && is.null(map))
    stop("Missing grid and map, please give at leat one.")
  
  if(!is.null(fixed_stations) && (length(intersect(class(fixed_stations),c("matrix","array","data.frame","SpatialPoints","SpatFD"))) == 0))
    stop("fixed_stations must be of class matrix, array, data.frame, SpatialPoints or SpatFD.")
  
  if(class(k) != "numeric" || length(k)==1){
    stop("k must be a positive integer.")
  }else if(round(k) != k){
    stop("k must be a positive integer.")
  }
  
  if(length(intersect(class(s0),c("matrix","array","data.frame","SpatialPoints"))) == 0){
    stop("s0 must be of class matrix, array, data.frame or SpatialPoints.")
  }
  
  if(ncol(s0)!=2)
    stop("s0 must have two columns.")
  
  if(length(intersect(class(model),c("variogramModel","list"))) == 0)
    stop("model must be a object of class variogramModel from gstat package or a lists of models of class variogramModel.")
  
  if(class(model)=="list"){
    if(!all(sapply(model,function(v)"variogramModel" %in% class(v))))
      stop("Every element in model must be of class variogramModel.")
  }
  
  if(!(method %in% c("lambda","scores")))
    stop("method must be one of 'lambda' or 'scores'.")
  
  if(!is.null(grid) && 
     length(intersect(c("matrix","array","data.frame"), class(grid))) == 0)
    stop("grid must be of class SpatialPoints, matrix, array or data.frame.")
  
  if(ncol(grid)!=2)
    stop("grid must have two columns.")
  
  if(is.null(grid) && !(inherits(map,"Spatial")))
    stop("map must be an Spatial object from sp package such as Line, Lines, Polygon, SpatialPolygons, SpatialGrid or SpatialPixels.")
  
  if("variogramModel" %in% class(model) && is.null(nharm)){
    message("As model is a single variogramModel and nharm is not specified, nharm is set to 1.")
    nharm = 1
  }
  
  
  # Basic objects
  if(is.null(grid))
    grid <- as.data.frame(sp::sample(map, n = 3e3,
                                     type = "regular"))
  
  G <- nrow(grid)
  
  stati0 <- grid[sample(1:G,k),]
  stati0 <- as.matrix(stati0)
  
  if("variogramModel" %in% class(model)){
    tmp_list <- list()
    for(i in 1:nharm)
      tmp_list[[i]] <- model
    model <- tmp_list
  }
  
  
  if(class(fixed_stations) == "SpatFD")
    fixed_stations <- fixed_stations$coords
  
  # Call optim
  stats::optim(par = c(stati0), fn = .vgm_model.fn,
               s0 = s0, model = model, method_i = method,
               grid = grid, fixed_stations = fixed_stations,
               method = "L-BFGS-B",
               control = list(
                 trace = 1L, 
                 factr = 1e-6,
                 REPORT = 4L),
               maxit = 100) -> result
  
  # Put into the grid the final result
  new_stati <- matrix(result$par,ncol=2)
  
  dist_tmp <- proxy::dist(new_stati,as.matrix(grid),diag = T)
  ids_tmp <- apply(as.matrix(dist_tmp),1,which.min)
  new_stati <- as.matrix(grid)[ids.tmp,]
  
  
  value <- list(new_stations = new_stati)
  
  # Plot -----------------------------------------
  if (plt){
    s0 <- as.data.frame(s0)
    
    if(is.null(map)){
      ggplot2::ggplot() +
        ggplot2::geom_point(ggplot2::aes(x = grid[,1],y = grid[,2]), colour = "gray90") + # grid
        ggplot2::geom_point(ggplot2::aes(x = s0[,1],y = s0[,2], col = "Target Points")) + # Target points
        ggplot2::geom_point(ggplot2::aes(x = new_stati[,1], y = new_stati[,2], col = " New Stations")) + # New stations
        if(!is.null(fixed_stations))
          ggplot2::geom_point(ggplot2::aes(x = fixed_stations[,1], y = fixed_stations[,2], col = "Fixed Stations")) + # Fixed stations
        ggplot2::theme_light() +
        ggplot2::labs(title = "Optimal Spatial Design",
                      subtitle = paste(k,"new stations")
                      x = "x", y = "y") + 
        ggplot2::theme(legend.position = ggplot2::element_blank()) -> final_plot
    } else {
      ggplot2::ggplot() +
        ggplot2::geom_point(ggplot2::aes(x = grid[,1],y = grid[,2]), colour = "gray90") + # grid
        ggplot2::geom_sf(data = sf::st_as_sf(map), fill = "ffffff00", colour = "gray40") + # map
        ggplot2::geom_point(ggplot2::aes(x = s0[,1],y = s0[,2], col = "Target Points")) + # Target points
        ggplot2::geom_point(ggplot2::aes(x = new_stati[,1], y = new_stati[,2], col = "New Stations")) + # New stations
        if(!is.null(fixed_stations))
          ggplot2::geom_point(ggplot2::aes(x = fixed_stations[,1], y = fixed_stations[,2], col = "Fixed Stations")) + # Fixed stations
        ggplot2::theme_light() +
        ggplot2::labs(title = "Optimal Spatial Design",
                      subtitle = paste(k,"new stations")
                      x = "Longitude", y = "Latitude") + 
        ggplot2::theme(legend.position = ggplot2::element_blank()) -> final_plot
    }
    
    value[["plot"]] <- final_plot
    
  }else{
    value[["plot"]] <- NULL
  }
  
  class(value) <- "OptimalSpatialDesign"
  return(value)
}


