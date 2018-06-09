#' Rao-Sampford sampling
#'
#' Draw a sample by means of Rao-Sampford sampling
#'
#' @param pik vector of first-order inclusion probabilities
#' @param n sample size
#' @param N population size (excluding self-selecting units)
#'
#' @note this function is a modified version of function \code{\link[sampling]{UPsampford}},
#' in the \pkg{sampling} package.
#'
#' @export

sampford <- function(pik, n, N, s, list){

    max_iter <- 5000
    sb = rep(2, N)
    y = pik/(1 - pik)/sum(pik/(1 - pik))
    step = 0
    while (sum(sb <= 1) != N & step <= max_iter) {
        sb <- as.vector(rmultinom(1, 1, pik/sum(pik)) +
                           rmultinom(1, .as_int(n - 1), y))
        step <- step + 1
    }
    if (sum(sb <= 1) == N){
        s[list] <- sb
        return(s)
    }else stop("Too many iterations. The algorithm was stopped.")
}
