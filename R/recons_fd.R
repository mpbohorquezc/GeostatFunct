
# Reconstruir objetos funcionales a partir de los scores y lambdas
recons_fd = function(X){

  if (inherits(X,"KS_pred")) {
    vari = X$name

    if(length(X) == 4) {

      if(inherits(X[[2]],"scores_pred")) {

        # Scores
        a = scores(X[[2]])
        mean_coef = X$SFD[[vari]]$fpca$meanfd$coefs
        nr=nrow(a)
        mean_coef = matrix(rep(mean_coef,nr),ncol = nr)

        coef_scores = (X$SFD[[vari]]$fpca$harmonics$coef %*% t(a)) + mean_coef
        result_scores = fd(coef_scores, X$SFD[[vari]]$fpca$harmonics$basis)

        result = result_scores

      } else {

        # Lambdas
        b = as.matrix(X[[2]]$lambda_pred)
        mean_coef = X$SFD[[vari]]$fpca$meanfd$coefs
        nr = ncol(b)
        mean_coef = matrix(rep(mean_coef,nr),ncol = nr)

        coef_lambda = (X$SFD[[vari]]$data_fd$coefs %*% b) + mean_coef
        result_lambda = fd(coef_lambda, X$SFD[[vari]]$data_fd$basis)

        result = result_lambda
      }

    } else {

      # Scores
      a = scores(X[[2]])
      mean_coef = X$SFD[[vari]]$fpca$meanfd$coefs
      nr=nrow(a)
      mean_coef = matrix(rep(mean_coef,nr),ncol = nr)

      coef_scores = (X$SFD[[vari]]$fpca$harmonics$coef %*% t(a)) + mean_coef
      result_scores = fd(coef_scores, X$SFD[[vari]]$fpca$harmonics$basis)

      # Lambdas
      b = as.matrix(X[[3]]$lambda_pred)

      mean_coef = X$SFD[[vari]]$fpca$meanfd$coefs
      nr=nrow(a)
      mean_coef = matrix(rep(mean_coef,nr),ncol = nr)

      coef_lambda = (X$SFD[[vari]]$data_fd$coefs %*% b) + mean_coef

      result_lambda = fd(coef_lambda, X$SFD[[vari]]$data_fd$basis)

      result = list(fd_scores = result_scores, fd_lambda = result_lambda)
    }



  }else{
    stop("Wrong class KS_pred")
  }
  return(result)
}
