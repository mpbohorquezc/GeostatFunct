print.SpatFD <- function(SpatFD, ...){

  data = SpatFD$data
  coords = SpatFD$coords
  data_fd = SpatFD$data_fd
  coefs = data_fd$coefs
  basis = data_fd$basis

  cat("## Data","\n")
  print(rbind(head(data, 3),
              rep("...", times=dim(data)[2]),
              tail(data, 3)))

  cat("\n","## Coords","\n")
  print(rbind(head(coords, 3),
              rep("...", times=dim(coords)[2]),
              tail(coords, 3)))

  cat("\n","## Data fd","\n")
  print(data_fd)

  cat("\n","## Coefs","\n")
  print(rbind(head(coefs, 3),
              rep("...", times=dim(coords)[2]),
              tail(coefs, 3)))

  cat("\n","## Basis","\n")
  print(head(basis))
}
}
