#' @title Differential Network Knockoff Filter
#'
#' @param AD AD sample data. n-p-q array or A list of length n, each of which is
#'     p-by-q spatial-temporal data matrix. n denotes the number of samples, p
#'     denotes the spatial dimension, q denotes the temporal dimension.
#' @param NC NC sample data. n-p-q array or A list of length n, each of which is
#'     p-by-q spatial-temporal data matrix. n denotes the number of samples, p
#'     denotes the spatial dimension, q denotes the temporal dimension.
#' @param screen_ratio Screen ratio. When screening, the number of screen is
#'     \eqn{ratio \cdot 0.5 \cdot p \cdot (p-1)}, (default: 0.5).
#' @param alpha Vector. The elasticnet mixing parameter, with \eqn{0\le\alpha\le 1}.
#' @param B Times of generating knockoff variables (default: 10)
#' @param fdr target false discovery rate (default: 0.1).
#' @param offset either 0 or 1 (default: 0).
#'
#' @return A list containing two lists, one of which is zero-one network and the
#'     other is the weight network
#' @import stats foreach
#' @export
#' @examples
#' # set.seed(2023)
#' # p<-50; q<-30; n1<-n2<-30
#' # AD <- vector(mode = "list", length = n1)
#' # NC <- vector(mode = "list", length = n2)
#' # for (i in seq(n1)) {
#' #   AD[[i]] <- matrix(rnorm(p*q), nrow = p, ncol = q, byrow = TRUE)
#' # }
#' # for (i in seq(n2)) {
#' #   NC[[i]] <- matrix(rnorm(p*q), nrow = p, ncol = q, byrow = TRUE)
#' # }
#' # result <- DNKF.filter(AD, NC, alpha = c(0.2, 0.3))
#'
#'
DNKF.filter <- function (AD, NC, alpha, screen_ratio = 0.5, B = 10, fdr = 0.1, offset = 0) {

  if (is.array(AD)) {
    AD_list <- vector(mode = "list", length = dim(AD)[1])
    for (i in seq(dim(AD)[1])) {
      AD_list[[i]] <- AD[i, , ]
    }
    AD <- AD_list
  }

  if (is.array(NC)) {
    NC_list <- vector(mode = "list", length = dim(NC)[1])
    for (i in seq(dim(NC)[1])) {
      NC_list[[i]] <- NC[i, , ]
    }
    NC <- NC_list
  }

  p <- dim(AD[[1]])[1]
  q <- dim(AD[[1]])[2]

  N1 <- length(AD)
  N2 <- length(NC)

  N_screen <- p*(p-1)/2 * screen_ratio

  cores <- parallel::detectCores()
  cl <- parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)

  AD_omega <- foreach(i = 1:N1, .errorhandling = 'pass', .packages = "DensParcorr") %dopar% {
    AD.clime <- DensParcorr::DensParcorr(t(AD[[i]]),dens.level =0.5,select=TRUE)
    AD_omega <- AD.clime$selected.precision
    rm(AD.clime)
    gc()
    return(AD_omega)
  }

  NC_omega <- foreach(i = 1:N2, .errorhandling = 'pass', .packages = "DensParcorr") %dopar% {
    NC.clime <- DensParcorr::DensParcorr(t(NC[[i]]),dens.level =0.5,select=TRUE)
    NC_omega <- NC.clime$selected.precision
    rm(NC.clime)
    gc()
    return(NC_omega)
  }

  X1_vec<-matrix(nrow=N1,ncol=p*(p-1)/2)
  W1_vec <- array()
  X2_vec<-matrix(nrow=N2,ncol=p*(p-1)/2)
  W2_vec <- array()

  W1<-lapply(lapply(AD_omega,cov2cor),mut_in_f)
  W2<-lapply(lapply(NC_omega,cov2cor),mut_in_f)

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
  X_screen <- getRho(X = X_vec, y = y, N_screen = N_screen)$X_screen
  screen_loc <- getRho(X = X_vec, y = y, N_screen = N_screen)$loc

  X <- X_screen
  network_01_list <- vector(mode = "list", length = length(alpha))
  network_weight_list <- vector(mode = "list", length = length(alpha))

  for (a in seq(length(alpha))) {
    foo <- knockoff::stat.lasso_coefdiff_bin
    statistic <- function(X, X_k, y) foo(X, X_k, y, alpha = alpha[a])
    knockoffs <- knockoff::create.second_order

    Xk <- foreach(i = seq(B), .errorhandling = 'stop', .packages = "knockoff") %dopar% {
      Xk <- knockoffs(X)
      return(Xk)
    }

    W <- foreach(i = seq(B), .combine = rbind, .errorhandling = 'stop', .packages = "knockoff") %dopar% {
      W <- statistic(X, Xk[[i]], y)
      return(W)
    }

    W_median <- c()
    for (i in seq(dim(W)[2])) {
      temp <- median(W[,i])
      W_median <- append(W_median, temp)
    }

    t <- knockoff::knockoff.threshold(W_median, fdr = fdr, offset = offset)

    selected <- sort(which(W_median >= t))
    threshold = t
    statistics <- W_median

    original_loc <- c()
    weight <- c()
    network_01 <- matrix(0, nrow = p, ncol = p)
    network_weight <- matrix(0, nrow = p, ncol = p)

    if (length(selected) != 0) {
      for (i in selected) {
        original_loc <- append(original_loc, screen_loc[i])
        weight <- append(weight, statistics[i])
      }
      original_loc <- getSubscripts(original_loc)
      for (i in seq(length(original_loc))) {
        network_01[original_loc[[i]][1], original_loc[[i]][2]] <- 1
        network_weight[original_loc[[i]][1], original_loc[[i]][2]] <- weight[i]
      }
    } else {
      cat(paste0("when alpha = ", alpha[a], ", select 0 edge"))
    }

    network_01_list[[a]] <- network_01
    network_weight_list[[a]] <- network_weight

  }

  parallel::stopCluster(cl)
  closeAllConnections()

  out <- list(network_01_list = network_01_list, network_weight_list = network_weight_list)
  return(out)
}
