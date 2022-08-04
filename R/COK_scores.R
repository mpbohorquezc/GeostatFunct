COK_scores <-
function(SFD,newcoords,model,vari=NULL,fill.all=T){

  #----------------------------------------------------------------------------
  #           VALIDANDO ARGUMENTOS *
  #----------------------------------------------------------------------------
  #all
  if(missing(SFD)){
    stop("Missing SFD")
  }
  if (missing(newcoords)){
    stop("Missing new coords")
  }
  if(missing(model)){
    stop("Missing model")
  }
  #SFD
  if(!inherits(SFD,"SpatFD")){
    stop("SFD must be an object SpatFD")
  }

  #newcoords
  if(!(is.matrix(newcoords) || is.data.frame(newcoords))){
    stop("Wrong class of newcoords object")
  }else if(!all(apply(newcoords, c(1,2), is.numeric))){
    stop("Newcoords must be numeric data")
  }else if(any(is.na(newcoords))){
    stop("There is some NA value in newcoords")
  }

  # messages default values
  if(missing(vari)){
    message("Using first variable by default")
  }
  if(missing(fill.all)){
    message("Using fill.all = TRUE by default")
  }

  #vari

  if(is.null(vari)){
    vari=1
  } else if ((is.character(vari)&& length(vari)==1)){
    if (length(which(names(SFD)==vari))==1){
      vari=which(names(SFD)==vari)
    }else if (length(which(names(SFD)==vari))==0){
      stop(paste(vari,"doesn't not exists. Change vari for an existing variable name."))
    }else if (length(which(names(SFD)==vari))==0){
      stop("There are more than one variable with the same name")
    }
  }
  if ((is.null(vari)  || !(is.numeric(vari)&& length(vari)==1))){
    stop("Wrong class of vari object")
  }

  #fill.all
  if ( !( ( isTRUE(fill.all) || isFALSE(fill.all) ) && length(fill.all)==1 ) ){
    stop("Wrong class of fill.all object")
  }


  #model

  if(!(inherits(model,"variogramModel") || inherits(model,"list"))){
    stop("Wrong class of model, model should be of class variogramModel or a list of them (use vgm of gstat package) ")
  }else if(inherits(model,"list") && !all(lapply(model,inherits,"variogramModel"))){
    stop("Wrong class of model, each element of list should be of class variogramModel (use vgm of gstat package)")
  }#else if(inherits(model,"list") && (length(model)!=ncol(as.data.frame(SFD[[vari]]$fpca$scores)))){
  #       stop("length of list of models must be equal to number of harmonics of the choosen variable ")
  #}else if(inherits(model,"variogramModel") && !(fill.all || (ncol(as.data.frame(SFD[[vari]]$fpca$scores))==1))){
  #        stop("If model is not a list and there are more than one nharm of that variable, then fill.all must be TRUE or you can create a list of models with the same number of harmonics")
  #}

  puntaje=list()
  puntajes=list()
  for(k in 1:length(SFD)){
    puntaje[[k]]=SFD[[k]]$fpca$scores
    rownames(puntaje[[k]])=SFD[[k]]$coordsname
    puntajes[[k]]=as.data.frame(puntaje[[k]])
    coordinates(puntajes[[k]])=SFD[[k]]$coords
  }

  colnames(newcoords)=c('x','y')
  coordinates(newcoords)=~x+y
  aa=rep(1:length(SFD),lapply(puntajes,ncol))
  bb=unlist(lapply(lapply(puntajes,ncol),seq))
  cc=paste0(
    c("g=gstat(,paste(colnames(SFD[[",rep("g=gstat(g,paste(colnames(SFD[[",(length(aa)-1))),
    aa,
    rep("]][[\"fpca\"]][[\"harmonics\"]][[\"coefs\"]])[",length(aa)),
    bb,
    rep("],",length(aa)),
    aa,
    rep(",sep=\".\"),puntajes[[",length(aa)),
    aa,
    rep("]][[",length(aa)),
    bb,
    rep("]]~1,puntajes[[",length(aa)),
    aa,
    rep("]])",length(aa))
  )

  eval(parse(text=cc))


  g <- gstat(g, model=model, fill.all=fill.all)

  vg <- variogram(g)
  mcl = fit.lmc(vg, g, fit.method=2)

  plot(vg, model = mcl)

  z = predict(mcl, newdata = newcoords)


  ret=list(SFD=SFD, COK_scores = z, model=mcl)
  class(ret)='COK_scores'
  return(ret)
}
