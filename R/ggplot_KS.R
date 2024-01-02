ggplot_KS <- function(
    KS, show.varpred = F, 
    main = "Functional Data", main2 = "Functional Data", 
    ylab = "Value", xlab = "Time", ndigits = 2, 
    palette.plot=c("#440154FF", "#3336FF", "#33FCFF", "#33FF4C", "#FDE725FF")){
    
    # Validation
    if(missing(KS)){
      stop("Missing KS")
    }
    
    if(!inherits(KS,"KS_pred")){
      stop("KS must be an object KS_pred")
    }
    
    # Functional object
    SFD <- SpatFD::recons_fd(KS)
    
    # Color palette
    custom_palette <- palette.plot
    
    # If just one method was used
    if(length(SFD) != 2){
      
      times <- SFD$basis$rangeval[1]:SFD$basis$rangeval[2]
      eval <- fda::eval.fd(times, SFD)
      melt_s <- suppressWarnings(reshape::melt(eval))
      melt_s$X2 <- as.factor(melt_s$X2)
      
      # Show predictions variance
      if(show.varpred){
        if(inherits(KS[[2]],"scores_pred")){
          levels(melt_s$X2) <- paste(levels(melt_s$X2), " (v_pred = ", round(KS[[2]]$scores_varpred$VTotal, ndigits), ")", sep = "")
          if(missing(main)){main <- "Functional Data - Scores Method"}
        } else {
          levels(melt_s$X2) <- paste(levels(melt_s$X2), " (v_pred = ", round(KS[[2]]$lambda_varpred$VTotal, ndigits), ")", sep = "")
          if(missing(main)){main <- "Functional Data - Lambda Method"}
        }
      }
      
      names(melt_s) = c("Time","Prediction","Value")
      
      # Color palette
      n_colors <- length(unique(melt_s$Prediction))
      color_palette <- colorRampPalette(custom_palette)(n_colors)
      
      # Plot
      graf=ggplot2::ggplot(melt_s, ggplot2::aes(x= Time, y= Value, col= Prediction)) +
        ggplot2::geom_line() +
        ggplot2::scale_color_manual(values = color_palette) +
        ggplot2::labs(title = main ) +
        ggplot2::labs(x = xlab, y = ylab) +
        ggplot2::theme_minimal()
      
      if(!show.varpred){
        graf <- graf + ggplot2::theme(legend.position="none")
      }
      
      return(graf)
      
    } else { # If both methods were used
      
      graf <- list()
      
      # Titles
      mainl <- list(main, main2)
      if(missing(main)){mainl[[1]] <- "Functional Data - Scores Method"}
      if(missing(main2)){mainl[[2]] <- "Functional Data - Lambda Method"}
      
      for(i in 1:2){
        
        times <- SFD[[i]]$basis$rangeval[1]:SFD[[i]]$basis$rangeval[2]
        eval <- fda::eval.fd(times, SFD[[i]])
        melt_s <- suppressWarnings(reshape::melt(eval))
        melt_s$X2 <- as.factor(melt_s$X2)
        
        # Show predictions variance
        if(show.varpred){
          if(inherits(KS[[1+i]],"scores_pred")){
            levels(melt_s$X2) <- paste(levels(melt_s$X2), " (v_pred = ", round(KS[[1+i]]$scores_varpred$VTotal, ndigits), ")", sep = "")
            
          } else {
            levels(melt_s$X2) <- paste(levels(melt_s$X2), " (v_pred = ", round(KS[[1+i]]$lambda_varpred$VTotal, ndigits), ")", sep = "")
            
          }
        }
        
        names(melt_s) = c("Time","Prediction","Value")
        
        # color palette
        n_colors <- length(unique(melt_s$Prediction))
        color_palette <- colorRampPalette(custom_palette)(n_colors)
        
        # Plot
        graf[[i]] =ggplot2::ggplot(melt_s, ggplot2::aes(x= Time, y= Value, col= Prediction)) +
          ggplot2::geom_line() +
          ggplot2::scale_color_manual(values = color_palette) +
          ggplot2::labs(title = mainl[[i]] ) +
          ggplot2::labs(x = xlab, y = ylab) +
          ggplot2::theme_minimal()
        
        if(!show.varpred){
          graf[[i]] <- graf[[i]] + ggplot2::theme(legend.position="none")
        }
        
      }
      
      return(graf)
    }
}
