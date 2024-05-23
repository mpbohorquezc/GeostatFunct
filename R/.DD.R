.DD <- function(expr, names, order = 1, debug=FALSE) {
  if (any(order>=1)) {  ## do we need to do any more work?
    w <- which(order>=1)[1]  ## find a derivative to compute
    if (debug) {
      cat(names,order,w,"\n")
    }
    ## update order
    order[w] <- order[w]-1
    ## recurse ...
    return(.DD(D(expr,names[w]), names, order, debug))
  }
  return(expr)
}
