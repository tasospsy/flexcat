estimateLCM <- function(X, ic = "aic"){
  # Estimating the latent class model
  f <- as.formula(paste0('cbind(',paste(names(X), collapse=','),')~1'))
  #f <- as.formula(paste("cbind(",paste(names(X)[1 : (J - 1)], "," , sep = "",collapse = ""), paste(names(X)[J], sep=""),") ~ 1", collapse=""))
  K <- 0L
  bestFit <- fit <- 1e100
  bestK <- 0
  model.fit <- list(AIC = rep(NA, 100), AICc = rep(NA, 100), AIC3 = rep(NA, 100), BIC = rep(NA, 100))

  X1 <- X - min(X) + 1L # if [0:..] transforms it to [1:...] for polca, else nothing happens and X1 = X
  repeat{
    K <- K + 1L
    cat('-> Estimating uLCM(',K,') at', format(Sys.time(), "%H:%M \n"))
    model.lc <- poLCA::poLCA(f, X1, nclass = K, maxiter = 5000)
    model.fit$AIC[K] <- model.lc$aic
    model.fit$AICc[K] <- model.lc$aic + (2 * model.lc$npar * (model.lc$npar + 1)) / (model.lc$Nobs - model.lc$npar - 1)
    model.fit$AIC3[K] <- -2 * model.lc$llik + 3 * model.lc$npar
    model.fit$BIC[K] <- model.lc$bic
    icNumber <- ifelse(sum(c("aic", "aicc", "aic3", "bic") %in% ic) != 1, 1, which(c("aic", "aicc", "aic3", "bic") %in% ic))
    #icNumber <- ifelse (ic %in% c("aic", "aicc", "aic3", "bic"), which(c("aic", "aicc", "aic3", "bic") == ic),1)
    fit <- model.fit[[icNumber]][K]
    if(fit < bestFit){
      bestFit <- fit
      bestK <- K
      bestModel <- model.lc
    }
    if (K == bestK + 2) break
  }

  ############################
  # Output of the latent class model
  # W are the class weights P(C = k)
  # P are the conditional response probabilities P(Xj = x|C = k): a K * nCat matrix for each item

  W <- bestModel$P
  P <- bestModel$probs
  return(list(W=W, P=P))
}
#estimateLCM(calib_set[1:J])
