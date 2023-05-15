#' @title Fisher Transformation
#'
#' @param x Data to be Fisher transformed.
#'
#' @return Data after Fisher transformation
#'
mut_in_f<-function(x){
  return(1/2*log((1+x)/(1-x)))
}
