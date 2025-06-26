
# LSCAT main administration function --------------------------------------
lscat_sim <- function(calibration_output,
                      precisionCriterion,
                      cutoff # only above this number, not =
) {
  resultsLcm <- calibration_output$resultsLcm
  A <- calibration_output$Sample_for_Sim
  X <- resultsLcm$X  # back to 0,1,.. format
  R <- resultsLcm$R  # back to 0,1,.. format
  RP <- resultsLcm$RP; RP <- RP/sum(RP)
  nR <- nrow(R)
  N <- nrow(X)
  M <- nrow(A)
  J <- ncol(X)
  nCat <- resultsLcm$nCat

  sr <- unique(rowSums(R)) # Andries wants SR out, but it helps here to define the classes
  YR <- g(R, intervals = list(min(sr):cutoff, (cutoff+1):max(sr)))
  YX <- g(X, intervals = list(min(sr):cutoff, (cutoff+1):max(sr)))
  YA <- g(A, intervals = list(min(sr):cutoff, (cutoff+1):max(sr)))

  Y <- sort(unique(YR))

  person <- 1L
  results <- list()
  cat("LSCAT in progress...\n")
  for (person in 1L : M){
    results[[person]] <- list()

    iter <- 0L
    obsPattern <- rep(FALSE, J)
    obsPatternScore <- matrix(NA, nR, J)

    PXi.ro  <- computePXi.ro(J, nCat, R, RP, obsPattern, obsPatternScore, iter)
    PY.ro   <- computePY.ro(J, R, RP, Y, YR, obsPattern, obsPatternScore, iter)
    PY.rox  <- computePY.rox(J, nCat, R, RP, Y, YR, obsPattern, obsPatternScore, iter)
    EY.ro   <- sum(PY.ro * Y)
    EY.rox  <- t(sapply(PY.rox, function(x) colSums(Y * x)))

    delta <- apply(PXi.ro * abs(EY.ro - EY.rox), 1, sum)
    selItem <- which.max(delta)[1]

    repeat{
      iter <- iter + 1L
      obsPattern[selItem] <- TRUE
      obsPatternScore[, selItem] <- A[person, selItem]
      PXi.ro  <- computePXi.ro(J, nCat, R, RP, obsPattern, obsPatternScore, iter)
      PY.ro   <- computePY.ro(J, R, RP, Y, YR, obsPattern, obsPatternScore, iter)

      if(iter == J | max(PY.ro) > precisionCriterion) {Yhat <- Y[which.max(PY.ro)]; break}

      PY.rox  <- computePY.rox(J, nCat, R, RP, Y, YR, obsPattern, obsPatternScore, iter)

      EY.ro   <- sum(PY.ro * Y)
      EY.rox  <- t(sapply(PY.rox, function(x) if (all(is.na(x))) return(rep(NA, nCat)) else return(colSums(Y * x))))
      delta <- apply(PXi.ro * abs(EY.ro - EY.rox), 1, sum, na.rm = TRUE)
      selItem <- which.max(delta)
    }

    #results[[person]]$observedItemScores   <- A[person, ]
    #results[[person]]$usedItemScores       <- obsPatternScore[1, ]
    results[[person]]$observedScore        <- YA[person]
    results[[person]]$estimatedScore       <- Yhat
    #results[[person]]$modalProbability     <- max(PY.ro)
    results[[person]]$efficiency           <- round((J - sum(obsPattern)) / J * 100, 1)
    #print(round(person / M * 100, 1)); Sys.sleep(.01); flush.console()
    progress <- 100*(person/M)
    if(progress %in% seq(10, 100, by = 10)) cat(progress,"%\n")
  }
  return(results)
}
