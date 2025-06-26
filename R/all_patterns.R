#' Generate All Response Patterns
#'
#' This function generates all possible response patterns for `J` items with `levels`.
#'
#' @param J Number of items
#' @param pat A vector for the item response categories
#'
#' @return A matrix of all possible response patterns (each row is one pattern)
#' @keywords internal
all.patterns <- function(J, pat){
  grid <- list()
  j <- 0;
  p <- length(pat)^J
  for (j in 1 : J){
    grid <- c(grid, j)
    grid[[j]] <- pat
  }
  X <- t(expand.grid(grid))
  dimnames(X) <- NULL
  return(X[J : 1, ])
}
