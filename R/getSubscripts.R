#' @title Get the location of the edges
#'
#' @param v vector denoting the location of the upper triangular matrix elements 
#'     after vectorization
#'
#' @return a list containing the location of each element in the matrix before 
#'     vectorization.
#'
getSubscripts <- function(v) {
  n <- 0
  res <- list()
  for (k in seq(length(v))) {
    while(plus(n) < v[k]) {
      n <- n + 1
    }
    j <- n + 1
    i <- v[k] - plus(n - 1)
    res[[k]] <- c(i, j)
  }
  return(res)
}