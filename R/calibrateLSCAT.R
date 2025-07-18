#' LSCAT Calibration
#' @description
#' This function first fits a Latent Class Model (LCM) using `poLCA::` and keeps the best fitting model according to a specified Information Criterion.
#' Then, computes the item-score density (π) under one of 4 choices: "mod", "emp", "inf" and "all"
#'
#' @param X A data frame with the available responses on all items
#' @param ic Information criterion; One of "aic", "aic3", "aicc", "bic". Default is "aic".
#' @param method Default is "all" (all possible responses).
#' Because the number of rows may be too large, several row-reduction methods have been added: "mod" (model-based responses), "emp" (empirical responses) and "inf" (most informative item responses).
#' @param saveLCMparams Logical. Default is `FALSE`. If `TRUE` the function returns also the class weights (W), and
#' the conditional response probabilities (P) under the best fitting Latent Class Model.
#'
#' @importFrom poLCA poLCA
#'
#' @references Psychogyiopoulos, A., Smits, N., Van der Ark, L.A.: Estimating the joint item- score density using an unrestricted latent class model: Advancing flexibility in computerized adaptive testing. Journal of Computerized Adaptive Testing (in press)
#' Van der Ark, L.A., Smits, N.: Computerized adaptive testing without IRT for flexible measurement and prediction. In: Essays on Contemporary Psychometrics, pp. 369–388. Springer, (2023). https://doi.org/10.1007/978-3-031-10370-4 19
#' @return An object of class list which includes
#' X : the data.frame used to calibrate the LCM.
#' R : A matrix with the item-score vectors used.
#' RP: A vector containing the probability of each item-score vector
#' nCat: The number of item categories
#' K: The number of classes in the best fitting LCM
#' @export
#'
#' @examples
calibrateLSCAT <- function(X, ic = "aic", method = "all", saveLCMparams = FALSE){

  # This function first fits a LCM using poLCA:: and keeps the best fitting model according to a specified IC
  # Then, computes the item-score density (π) under one of 3 choices: "mod", "emp", and "all)
  ############################
  # Preliminaries
  # J = number of items
  # ic = information criterion: "aic", "aic3", "aicc", "bic"
  # nCat = number of categories
  # Checks of class, missing values and item scores

  J <- ncol(X)
  nCat <- max(X) - min(X) + 1L

  if(!is.data.frame(X)) stop('X is not a data.frame')
  if(any(is.na(X))) stop('Some item scores are missing')
  #if(!all(unlist(X) %in% 1:2)) stop('Data contain item scores other than 1 and 2')

  # -------------------------------------------------------------------------
  .lcmout <- estimateLCM(X, ic) # see estimateLCM.R
  W <- .lcmout$W
  P <- .lcmout$P
  bestK <- length(W)
  # -------------------------------------------------------------------------

  ############################
  # p.method
  # Managing the output of the latent class model
  # nR <- nCat^J is the number of possible response patters. If J is large, nR may become infeasibly large
  # R is an nR x J matrix of all response patterns
  # Because the number of rows may be too large, several row-reduction methods  have been added: "emp" (empirical), "mod" (model based)
  # RP is an nR x 1 vector of the probability that a response patterns occurs under the latent class model
  # Note: R is based on X1, the dataset used for poLCA
  X1 <- X - min(X) + 1L
  if (method == "all") R <- t(all.patterns(J, 1 : nCat))

  if (method == "emp") {
    Rchar <- matrix(sort(unique(apply(X1, 1, paste0, collapse = ""))))
    R <- t(sapply(Rchar[,1], function(s) as.numeric(strsplit(s, "")[[1]])))
    # @tofix 'Rchar' has always min val 1
  }
  if (method == "mod") {
    Nsim <- 10000                                                                                     # Number of simulees
    Ksim <- sample(1:bestK, Nsim, TRUE, prob= W)                                                      # Classes of the simulees
    Psim <- lapply(Ksim, function(n) do.call(rbind, lapply(P, function(mat) mat[n, , drop = FALSE]))) # Response probablities of the simulees
    Psim <- lapply(Psim, function(x) t(apply(x, 1, cumsum)))                                          # Cumulative response probablities of the simulees
    Usim <- lapply(1:Nsim, function(r) matrix(runif(J), nrow = J, ncol = nCat, byrow = FALSE))        # Random numbers
    PgtU <- mapply(function(p, u) (p >= u) * 1, Psim, Usim, SIMPLIFY = FALSE)                         # Compare Cum resp probabilities to random numbers
    Xsim <- matrix(unlist(lapply(PgtU, function(x) apply(x, 1, sum))), Nsim, J, byrow = TRUE)         # Simulated item score data
    Rchar <- matrix(sort(unique(apply(Xsim, 1, paste0, collapse = ""))))                                 # see "emp"
    R <- t(sapply(Rchar[,1], function(s) as.numeric(strsplit(s, "")[[1]])))                           # see "emp"
  }

  # -------------------------------------------------------------------------
  # Compute π
  nR <- nrow(R)
  RP <- rep(NA, nR)
  for (p in 1 : nR){
    RP[p] <- 0
    for (k in 1 : bestK){
      RP.k <- log(W[k])
      for(j in 1 : J) RP.k <- RP.k + log(P[[j]][k, R[p, j]])
      #cat(P[[j]][k, R[p, j]])
      RP[p] <- RP[p] + exp(RP.k)
    }
  }
  RP <- RP/sum(RP) # normalize π
  # return the data to original format
  if(min(X1)-min(X)==1) R <- R - 1L

  result <- list(X = X, R = R, RP = RP, nCat = nCat, K = bestK)

  if (saveLCMparams) {
    result$W <- W
    result$P <- P
  }
  return(result)
}
