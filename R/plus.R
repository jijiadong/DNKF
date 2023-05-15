#' @title Accumulation function
#'
#' @param n a positive integer.
#'
#' @return a positive integer is equal to the accumulation from 1 to n.
#'
plus <- function(n) {
  if (n == 0) {
    return(1)
  } else return(n + plus(n-1))
}
