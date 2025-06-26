#' Compute the Item-score Density
#' @description
#' This function computes the item-score density (π), i.e., the vector including the probabilities
#' of each score ppatern. It uses the output of the estimated LCM from the function
#' flexcat::estimateLCM()
#'
#' @param R A matrix of possible response patterns
#' @param W Class weight probabilities: an ULCM parameter
#' @param P Item-response probabilities given class membership: an ULCM parameter
#' @param bestK The number of the latent classes that best describe ULCM fit.
#' @param J The number items 'J'
#'
#' @references Psychogyiopoulos, A., Smits, N., Van der Ark, L.A.: Estimating the joint item- score density using an unrestricted latent class model: Advancing flexibility in computerized adaptive testing. Journal of Computerized Adaptive Testing (in press)
#' Van der Ark, L.A., Smits, N.: Computerized adaptive testing without IRT for flexible measurement and prediction. In: Essays on Contemporary Psychometrics, pp. 369–388. Springer, (2023). https://doi.org/10.1007/978-3-031-10370-4 19
#' @return A vector with probabilities that sum to 1.

#' @export
#'
#' @examples
#'
compute_p <- function(R, bestK=NULL, W, P, J=NULL){
## It needs R in 1:2 format!
if(min(R)==0) R <- R + 1
if(is.null(bestK)) bestK <- length(W)
if(is.null(J)) J <- ncol(R)
nR <- nrow(R)
RP <- rep(NA, nR)
for (p in 1 : nR){
  RP[p] <- 0
  for (k in 1 : bestK){
    RP.k <- log(W[k])
    for(j in 1 : J) RP.k <- RP.k + log(P[[j]][k, R[p, j]])
    RP[p] <- RP[p] + exp(RP.k)
  }
}
RP <- RP/sum(RP) # normalize π
return(RP)
}
