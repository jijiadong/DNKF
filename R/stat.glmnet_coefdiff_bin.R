#' @title Importance statistics based on a GLM with cross-validation
#'
#' @param X n-by-p matrix of original variables.
#' @param X_k n-by-p matrix of knockoff variables.
#' @param y vector of length n, containing the response variables.
#' @param Zv a vector or matrix, containing the covariates.
#' @param family response type.The default response family is 'binomial', for a
#'  logistic regression model.
#' @param ... Other arguments that can be passed to glmnet.
#' @return A vector of statistics \eqn{W} of length p.
#'
stat.glmnet_coefdiff_bin <- function(X, X_k, Zv, y, family='binomial',...) {

  # Randomly swap columns of X and Xk
  swap <- rbinom(ncol(X),1,0.5)
  swap.M <- matrix(swap,nrow=nrow(X),ncol=length(swap),byrow=TRUE)
  X.swap  <- X * (1-swap.M) + X_k * swap.M
  Xk.swap <- X * swap.M + X_k * (1-swap.M)

  p <- ncol(X)

  cv_coeffs_glmnet <- function(X, y, family = 'binomial', nlambda=500, intercept=T, parallel=F,...) {
    # Standardize variables
    X <- scale(X)

    n <- nrow(X); p <- ncol(X)

    if (!methods::hasArg(family) ) family <- 'binomial'
    if (!methods::hasArg(lambda) ) {
      if( identical(family, "gaussian") ) {
        if(!is.numeric(y)) {
          stop('Input y must be numeric.')
        }
        # Unless a lambda sequence is provided by the user, generate it
        lambda_max <- max(abs(t(X) %*% y)) / n
        lambda_min <- lambda_max / 2e3
        k <- (0:(nlambda-1)) / nlambda
        lambda <- lambda_max * (lambda_min/lambda_max)^k
      }
      else {
        lambda <- NULL
      }
    }

    alpha_seq <- seq(0, 1, 0.04)

    cv.glmnet.all <- glmnetUtils::cva.glmnet(X, y, family=family, lambda=lambda, alpha = alpha_seq, intercept=intercept,
                                             standardize=F,standardize.response=F, parallel=parallel,...)

    get_best <- function(model_all) {
      get_cvm <- function(model) {
        index <- match(model$lambda.min, model$lambda)
        model$cvm[index]
      }
      enet_performance <- data.frame(alpha = model_all$alpha)
      models <- model_all$modlist
      enet_performance$cvm <- vapply(models, get_cvm, numeric(1))
      minix <- which.min(enet_performance$cvm)
      best_alpha <- model_all$alpha[minix]
      return(best_alpha)
    }

    best_alpha <- get_best(cv.glmnet.all)
    cv.glmnet.fit <- glmnet::cv.glmnet(X, y, family=family, lambda=lambda, alpha=best_alpha, intercept=intercept,
                                       standardize=F,standardize.response=F, parallel=parallel,...)
    coef(cv.glmnet.fit, s = "lambda.min")
  }

  # Compute statistics
  glmnet.coefs <- cv_coeffs_glmnet(cbind(X.swap, Xk.swap, Zv), y, family=family,...)
  if(family=="multinomial") {
    Z <- abs(glmnet.coefs[[1]][2:(2*p+1)])
    for(b in 2:length(glmnet.coefs)) {
      Z <- Z + abs(glmnet.coefs[[b]][2:(2*p+1)])
    }
  } else if (family=="cox") {
    Z <- glmnet.coefs[1:(2*p)]
  } else {
    Z <- glmnet.coefs[2:(2*p+1)]
  }
  orig <- 1:p
  W <- abs(Z[orig]) - abs(Z[orig+p])

  # Correct for swapping of columns of X and Xk
  W <- W * (1-2*swap)

  return(W)
}

