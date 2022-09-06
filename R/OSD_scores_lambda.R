.vgm_model.fn <- function(points, s0, model, method_i = "lambda", grid){
  
  # validation ---------------------------------
  
  # This is an internal function, thus its arguments will be validated in the
  # main function FD_optimal_design.
  
  # Basic objects
  method <- method_i
  N <- nrow(s0)
  k <- length(points)/2
  z <- rep.int(1,k)
  L <- length(model)
  points <- matrix(points,ncol=2)
  colnames(points) <- c("x","y")
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
    
    tot_variance <- 0
    for(l in 1:N){
      #c_vec <- rep(0,k)
      for(j in 1:L){
        #c_vec <- c_vec + gstat::variogramLine(model[[j]], dist_vector = dist_s0[,l],coveriance = T)
        c_vec <- gstat::variogramLine(model[[j]], dist_vector = dist_s0[,l],coveriance = T)
        tot_variance <- c_vec%*%omegas[[j]]%*%t(c_vec)
        # Sugerencia: Ponderar por proporcion de varianza explicada.
      }
    }
    
    return(tot_variance)
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

FD_optimal_design <- function(k,s0,model,method = "lambda", grid = NULL,
                              map = NULL,nharm = NULL){
  # validation ----------------------------------------------
  
  if (missing(k))
    stop("Missing k.")
  
  if(missing(s0))
    stop("Missing s0.")
  
  if(missing(s0))
    stop("Missing Model.")
  
  if(is.null(grid) && is.null(map))
    stop("Missing grid and map, please give at leat one.")
  
  
  
  
  if(class(k) != "numeric" || length(k)==1){
    stop("k must be a positive integer.")
  }else if(round(k) != k){
    stop("k must be a positive integer.")
  }
  
  if(length(intersect(class(s0),c("matrix","array","data.frame","SpatialPoints"))) == 0){
    stop("s0 must have class matrix, array, data.frame or SpatialPoints.")
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
  
  # Call optim
  stats::optim(par = c(stati0), fn = .vgm_model.fn,
               s0 = s0, model = model, method_i = method,
               grid = grid,
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
  
  return(new_stati)
}

