#' @title Calculate Pearson Correlation Coefficient of X and y
#'
#' @param X n-by-p matrix of predictors.
#' @param y response vector of length n.
#' @param N_screen the number of columns left.
#'
#' @return a list containing a location vector denoting the location of the 
#'     screened out columns and the Matrix after screening.
#'
getRho <- function (X, y, N_screen) {
  p <- dim(X)[2]
  rho <- c()
  for (i in seq(p)) {
    temp <- cor(X[, i], y)
    rho <- append(rho, temp)
  }
  rho <- abs(rho)
  threshold <- sort(rho, decreasing = TRUE)[N_screen]
  loc <- which(rho >= threshold)
  X_screen <- X[, loc]
  
  out <- list(loc = loc, X_screen = X_screen)
  return(out)
}
