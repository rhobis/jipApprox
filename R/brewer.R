#' Brewer sampling procedure --------------------------------------------------
#'
#' @param pik vector of first-order inclusion probabilities
#' @param n sample size
#' @param N population size
#' @param s vector of length N, with 1s at the positions of self-selecting units
#' @param list vector with positions of self selcting units
#'
#'
#' @note this function is a modified version of function \code{\link[sampling]{UPbrewer}},
#' from the \pkg{sampling} package.
#'
#' @export

brewer <- function(pik, n, N, s, list){
    sb <- rep(0, N)
    for (j in 1:n) {
        a <- sum(pik * sb)
        p <- (1 - sb) * pik * ((n - a) - pik)/((n - a) - pik * (n - j + 1))
        p <- p/sum(p)
        p <- cumsum(p)
        u <- runif(1)
        # for (h in seq_along(p)){
        #     if (u < p[h]) break
        # }
        h <- which(u<p)[1]
        sb[h] <- 1
    }
    s[list] <- sb
    return(s)
}
