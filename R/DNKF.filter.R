#' @title Differential Network Knockoff Filter
#'
#' @param screen_ratio Screen ratio. When screening, the number of screen is
#'     \eqn{ratio \cdot 0.5 \cdot p \cdot (p-1)}, (default: 0.5).
#' @param B Times of generating knockoff variables (default: 10)
#' @param fdr target false discovery rate (default: 0.1).
#' @param offset either 0 or 1 (default: 0).
#' @param data1 One group data. p-q-n1 array or A list of length n1, each of which is
#'     p-by-q spatial-temporal data matrix. n1 denotes the number of samples, p
#'     denotes the spatial dimension, q denotes the temporal dimension.
#' @param data2 Another group data. p-q-n2 array or A list of length n2, each of which is
#'     p-by-q spatial-temporal data matrix. n2 denotes the number of samples, p
#'     denotes the spatial dimension, q denotes the temporal dimension.
#' @param Zv1 A n1*v covariates matrix corresponding to data1.
#' @param Zv2 A n2*v covariates matrix corresponding to data2.
#' @param cores Number of cores used in parallel (default: Number of all cores
#'     of the computer).
#'
#' @return Differential network matrix. The non-zero elements of the matrix represent
#'     the estimated difference of partial correlation coefficients
#' @import stats foreach
#' @export
#' @examples
#' p <- 6; q <- 20; n <- 15
#' set.seed(2023)
#' Theta <- diag(p)
#' Theta[1,2] <- Theta[1,3] <- Theta[4,5] <- Theta[4,6] <- 1
#' Theta <- t(Theta)+Theta
#' diag(Theta) <- 1
#' omega1 <- Theta * sample(c(-1, 1), p * p, replace = TRUE) * runif(p * p, 0.4, 0.5)
#' omega1[lower.tri(omega1, diag = TRUE)] <- 0
#' omega1 <- as.matrix(omega1)
#' omega1 <- omega1 + t(omega1)
#' diag(omega1) <- abs(min(eigen(omega1)$values)) + 1
#' sigma1 <- cov2cor(solve(omega1))
#' omega1 <- solve(sigma1)
#' omega1[abs(omega1)<10^-4] <- 0
#' omega2 <- omega1
#' omega2[1:(p/2),1:(p/2)] <- -1*omega2[1:(p/2),1:(p/2)]
#' diag(omega2) <- diag(omega1)
#' sigma2 <- solve(omega2)
#' delta <- cov2cor(omega1)-cov2cor(omega2)
#' sigmaT1 <- 0.4^abs(outer(1:q,1:q,"-"))
#' sigmaT2 <- 0.5^abs(outer(1:q,1:q,"-"))
#'
#' rmatnorm <- function(n,mean,sigmaS,sigmaT){
#'   nr <- nrow(sigmaS)
#'   nc <- ncol(sigmaT)
#'   R <- chol(sigmaT, pivot = TRUE)
#'   R <- R[, order(attr(R, "pivot"))]
#'   Q <- chol(sigmaS, pivot = TRUE)
#'   Q <- Q[, order(attr(Q, "pivot"))]
#'   mat = array(dim = c(nr, nc, n))
#'   for (i in 1:n) {
#' 	   mat[, , i] = mean + t(Q) %*% matrix(rnorm(nr * nc), nrow = nr) %*% R
#'   }
#'   mat
#' }
#'
#' X1 <- rmatnorm(n,mean=matrix(0,p,q),sigmaS=sigma1,sigmaT=sigmaT1)
#' X2 <- rmatnorm(n,mean=matrix(0,p,q),sigmaS=sigma2,sigmaT=sigmaT2)
#'
#' result <- DNKF.filter(data1 = X1, data2 = X2,  B = 4, cores = 2)
#' result
#'
DNKF.filter <- function (data1, data2, Zv1=NULL, Zv2=NULL, screen_ratio = 0.5, B = 10, fdr = 0.1, offset = 0, cores) {

  if (is.array(data1)) {
    data1_list <- vector(mode = "list", length = dim(data1)[3])
    for (i in seq(dim(data1)[3])) {
      data1_list[[i]] <- data1[, ,i ]
    }
    data1 <- data1_list
  }

  if (is.array(data2)) {
    data2_list <- vector(mode = "list", length = dim(data2)[3])
    for (i in seq(dim(data2)[3])) {
      data2_list[[i]] <- data2[, ,i ]
    }
    data2 <- data2_list
  }

  p <- dim(data1[[1]])[1]
  q <- dim(data1[[1]])[2]

  N1 <- length(data1)
  N2 <- length(data2)

  N_screen <- p*(p-1)/2 * screen_ratio
  N_screen <- floor(N_screen)

  if (!is.null(Zv1) & !is.null(Zv2)) {
    if (is.vector(Zv1) & is.vector(Zv2)){
	  if (length(Zv1)!=N1 | length(Zv2)!=N2)
        stop("The length of Zv1 or Zv2 is incorrect!")
	  Zv <- as.matrix(c(Zv1,Zv2))
	} else {
	  if (nrow(Zv1)!=N1 | nrow(Zv2)!=N2)
        stop("The dimension of Zv1 or Zv2 is incorrect!")
      Zv <- rbind(Zv1,Zv2)
	  Zv <- as.matrix(Zv)
	}
  } else{
    Zv <- NULL
  }


  if (missing(cores)) {
    cores <- parallel::detectCores()
  }
  cl <- parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)

  data1_omega <- foreach(i = 1:N1, .errorhandling = 'stop', .export=c("DensPcorr"),.packages = "DensParcorr") %dopar% {
    data1.clime <- DensPcorr(t(data1[[i]]),dens.level =0.5,select=TRUE)
    data1_omega <- data1.clime$selected.precision
    rm(data1.clime)
    gc()
    return(data1_omega)
  }

  data2_omega <- foreach(i = 1:N2, .errorhandling = 'stop',.export=c("DensPcorr"), .packages = "DensParcorr") %dopar% {
    data2.clime <- DensPcorr(t(data2[[i]]),dens.level =0.5,select=TRUE)
    data2_omega <- data2.clime$selected.precision
    rm(data2.clime)
    gc()
    return(data2_omega)
  }

  X1_vec <- matrix(nrow=N1,ncol=p*(p-1)/2)
  X2_vec <- matrix(nrow=N2,ncol=p*(p-1)/2)

  mut_in_f<-function(x){
    return(1/2*log((1+x)/(1-x)))
  }

  R1 <- lapply(data1_omega,cov2cor)
  R2 <- lapply(data2_omega,cov2cor)

  R1_average <- matrix(0, nrow = dim(R1[[1]])[1], ncol = dim(R1[[1]])[2])
  R2_average <- matrix(0, nrow = dim(R2[[1]])[1], ncol = dim(R2[[1]])[2])

  for (mat in R1) {
    R1_average <- R1_average + mat
  }
  R1_average <- R1_average / length(R1)
  for (mat in R2) {
    R2_average <- R2_average + mat
  }
  R2_average <- R2_average / length(R2)

  R_delta_average <- R1_average - R2_average

  W1<-lapply(R1, mut_in_f)
  W2<-lapply(R2, mut_in_f)

  #Generate data for logistic regression
  for (i in 1:N1) {
    W1_vec <- as.vector(W1[[i]][upper.tri(W1[[i]], diag = FALSE)])
    X1_vec[i,] <- W1_vec
  }
  for (i in 1:N2) {
    W2_vec <- as.vector(W2[[i]][upper.tri(W2[[i]], diag = FALSE)])
    X2_vec[i,] <- W2_vec
  }

  X_vec <- rbind(X1_vec, X2_vec)
  y1 <- rep(1, N1)
  y2 <- rep(0, N2)
  y <- append(y1, y2)

  # Screen
  getRho <- function (X, y, N_screen) {

	rho <- abs(cor(X, y))
    threshold <- sort(rho, decreasing = TRUE)[N_screen]
    loc <- which(rho >= threshold)
    X_screen <- X[, loc]

    out <- list(loc = loc, X_screen = X_screen)
    return(out)
  }

  X_screen <- getRho(X = X_vec, y = y, N_screen = N_screen)$X_screen
  screen_loc <- getRho(X = X_vec, y = y, N_screen = N_screen)$loc

  X <- X_screen
  n <- nrow(X)
  stopifnot(length(y) == n)

  W_list <- foreach(i = seq(B), .errorhandling = 'stop', .export=c("stat.glmnet_coefdiff_bin"),
                    .packages = "knockoff") %dopar% {
    Xk <- knockoff::create.second_order(X)
    W_list <- stat.glmnet_coefdiff_bin(X, Xk, Zv, y)
    gc()
    return(W_list)
  }

  parallel::stopCluster(cl)
  closeAllConnections()

  W_mat <- do.call(rbind, W_list)

  W_median <- robustbase::colMedians(W_mat)

  t_median <- knockoff::knockoff.threshold(W_median, fdr=fdr, offset=offset)
  selected_median <- sort(which(W_median >= t_median))

  threshold = t_median
  statistics <- W_median
  selected <- selected_median

  original_loc <- c()
  network_01 <- matrix(0, nrow = dim(R_delta_average)[1], ncol = dim(R_delta_average)[2])

  matrix_temp <- matrix(0, nrow = dim(R_delta_average)[1], ncol = dim(R_delta_average)[2])
  indices <- which(upper.tri(matrix_temp), arr.ind = TRUE)
  rm(matrix_temp)

  if (length(selected) != 0) {
    for (i in selected) {
      original_loc <- append(original_loc, screen_loc[i])
    }
    original_loc_two <- indices[original_loc, ]
    network_01[original_loc_two] <- 1
    network_01 <- network_01 + t(network_01)
  } else {
    cat("0 edges are selected")
  }
  network_weight <- network_01 * R_delta_average
  return(network_weight)
}
