#' @title Dens-based approach for precision matrix estimation.
#'
#' @details
#' This function implements the statistical method proposed in Wang et al. (2016).
#'     See Rpackage "DensParcorr".
#'
#' @param data Input data matrix of size n (observations) times p (variables).
#' @param select Whether to conduct the Dens-based selection. If FALSE, output
#'     will only contain the estimated partial correlation list and precision
#'     matrix list corresponding to the default tuning parameter series ranging
#'     from 1e-8 to 0.6. If TRUE, the ouput will include the previous results and
#'     the selected partial correlation matrix and percision matrix corresponding
#'     to the specified density level. Default is FALSE.
#' @param dens.level Specify the density level in Dens-based tuning parameter selection
#'     method (0<dens.level<1). This option is valid only when select=TRUE.
#' @param plateau.thresh The criterion to select the plateau.
#' @param Parcorr.est Previous output from DensPcorr function.
#' @param lambda The tuning parameters for estimating the precision matrix
#'     ranging from 0 to 1.
#'
#' @return An R list containing the following terms: \item{selected.precision}{Selected Precision matrix corresponding to dens.level.} \item{selected.lambda}{Selected tuning parameter corresponding to dens.level.}
DensPcorr <- function (data, select = FALSE, dens.level = 0.5,
                       plateau.thresh = 0.01, Parcorr.est = NULL, lambda = NULL) {

  if (!is.null(lambda)) {
    Prec.mat = clime::clime(data, lambda = lambda)$Omegalist[[1]]
    Par.mat = DensParcorr::prec2part(Prec.mat)

    Results = list()
    Results$selected.partial.corr = Par.mat
    Results$selected.precision = Prec.mat
    Results$selected.lambda = lambda
    Results$selection.method = paste("Lambda =", lambda)
    return(Results)
  }
  lambda.max = 0.6
  lambda.min = 1e-08
  nlambda = 10
  lambda = 10^(seq(log10(lambda.min), log10(lambda.max), length.out = nlambda))
  data = as.matrix(data)
  if (!is.null(Parcorr.est)) {
    if (is.null(Parcorr.est$precision.list) | is.null(Parcorr.est$lambda))
      stop("Parcorr.est is modified!")
    Prec.mat = Parcorr.est$precision.list
    lambda = Parcorr.est$lambda
    Prec.mat = Prec.mat[order(lambda)]
    lambda = lambda[order(lambda)]
  } else {
    lambda = lambda[order(lambda)]
    Prec.mat = clime::clime(data, lambda = lambda)$Omegalist
  }
  dens = vector()
  for (i in 1:length(Prec.mat)) {
    dens[i] = DensParcorr::prec2dens(Prec.mat[[i]])
  }
  while (abs(1 - sort(dens)[length(dens) - 1]/max(dens)) > 0.05) {
    lambda = c(min(lambda)/10, lambda)
    Prec.mat = append(clime::clime(data, lambda = lambda[1])$Omegalist,
                      Prec.mat)
    dens = c(DensParcorr::prec2dens(Prec.mat[[1]]), dens)
  }
  if (select) {
    if (dens.level == "plateau") {
      select.index = max(which((1 - dens/max(dens)) <=
                                 plateau.thresh))
      while (dens[select.index] == max(dens)) {
        lam.min = lambda[select.index]
        lam.max = lambda[select.index + 1]
        lambda2 = 10^(seq(log10(lam.min), log10(lam.max),
                          length.out = 4))[2:3]
        Prec.mat2 = clime::clime(data, lambda = lambda2)$Omegalist
        Prec.mat = append(Prec.mat2, Prec.mat)
        lambda = c(lambda2, lambda)
        Prec.mat = Prec.mat[order(lambda)]
        lambda = lambda[order(lambda)]
        for (i in 1:length(Prec.mat)) {
          dens[i] = DensParcorr::prec2dens(Prec.mat[[i]])
        }
        select.index = max(which((1 - dens/max(dens)) <= plateau.thresh))
      }
    } else if (is.numeric(dens.level) & dens.level < 1 & dens.level > 0) {
      select.index = which.min(abs(dens - max(dens) * dens.level))
      while (abs(dens.level - dens[select.index]/max(dens)) > 0.05) {
        if (max(dens) * dens.level > dens[select.index]) {
          lam.min = lambda[select.index - 1]
          lam.max = lambda[select.index]
        } else {
          lam.min = lambda[select.index]
          lam.max = lambda[select.index + 1]
        }
        lambda2 = 10^(seq(log10(lam.min), log10(lam.max),
                          length.out = 4))[2:3]
        Prec.mat2 = clime::clime(data, lambda = lambda2)$Omegalist
        Prec.mat = append(Prec.mat2, Prec.mat)
        lambda = c(lambda2, lambda)
        Prec.mat = Prec.mat[order(lambda)]
        lambda = lambda[order(lambda)]
        for (i in 1:length(Prec.mat)) {
          dens[i] = DensParcorr::prec2dens(Prec.mat[[i]])
        }
        select.index = which.min(abs(dens - max(dens) * dens.level))
      }
    }
  }

  Par.mat = list()
  for (i in 1:length(lambda)) {
    Par.mat[[i]] = DensParcorr::prec2part(Prec.mat[[i]])
  }
  Results = list()
  if (select) {
    Results$selected.partial.corr = Par.mat[[select.index]]
    Results$selected.precision = Prec.mat[[select.index]]
    Results$selected.lambda = lambda[[select.index]]
    Results$partial.corr.list = Par.mat
    Results$precision.list = Prec.mat
    Results$lambda.list = lambda
  }
  else {
    Results$partial.corr.list = Par.mat
    Results$precision.list = Prec.mat
    Results$lambda.list = lambda
    Results$Dens = dens
    Results$Dens.Percentage = round(dens/max(dens), 3)
  }
  if (select) {
    if (dens.level == "plateau") {
      Results$selection.method = "Dens-plateau"
    } else if (is.numeric(dens.level) & dens.level < 1 & dens.level > 0) {
      Results$selection.method = paste(round(dens.level *
                                               100, 1), "% Dens (Actual=", round(dens[select.index]/max(dens) *
                                                                                   100, 1), "%)", sep = "")
    }
    Results$Dens = dens
    Results$Dens.Percentage = round(dens/max(dens), 3)
  }
  return(Results)
}
