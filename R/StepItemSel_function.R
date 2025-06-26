.StepItemSel <- function(X){
  steps <- list()
  # Step 1
  step <- 1
  init.cases <- data.frame()
  for(i in 1:J) {
    for(j in 1:J) {
      if(i != j) init.cases <- rbind(init.cases, data.frame(Y = i, X1 = j))
    }
  }
  plot(0, 0, type = "n", xlim = c(1, nrow(init.cases)), ylim = c(0, 2000),
       xlab = "case", ylab = "AIC")
  aic <- c()
  for(row in 1:nrow(init.cases)){
    f <- as.formula(paste0(names(X)[init.cases[row,1]],' ~ ',names(X)[init.cases[row,2]]))
    model <- glm(f, family = binomial, data = X)
    points(row, AIC(model))
    aic <- c(aic,AIC(model))
  }
  steps[[step]] <- list(AIC=min(aic), vars = init.cases[which.min(aic),])
  steps[[step]]$formula <- paste0(names(X)[steps[[step]]$vars[,'Y']],' ~ ',
                                  names(X)[steps[[step]]$vars[,paste0("X",step)]])
  # Step 2+
  # Now we add one more item per step and check model fit
  repeat({
    step <- step + 1
    available_vars <- if (step == 2) {
      tmp <- 1:J
      tmp[-unlist(steps[[step - 1]]$vars)]
    } else {
      available_vars[available_vars != steps[[step - 1]]$vars]
    }
    iter <- 0; aic <- c()
    for (var in available_vars) {
      iter <- iter + 1
      f <- as.formula(paste0(steps[[step - 1]]$formula, " + ", names(X)[var]))
      model <- glm(f, family = binomial, data = X)
      points(iter, AIC(model))
      aic <- c(aic, AIC(model))
    }
    steps[[step]] <- list(AIC = min(aic),vars = available_vars[which.min(aic)]
    )
    steps[[step]]$formula <- paste0(steps[[step - 1]]$formula," + ",
                                    names(X)[steps[[step]]$vars])
    if (steps[[step - 1]]$AIC <= steps[[step]]$AIC) break
  })
  best <- steps[[step-1]]$formula
  items <- unlist(strsplit(gsub("~|\\+", " ", best), "\\s+"))
  #item_pattern <- names(X[1:J]) %in% items
  return(items)
  dev.off()
}
