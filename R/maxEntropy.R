#' Conditional Poisson Sampling  - compute selection probabilities
#'
#' Compute matrix of selection probabilities for Conditional Poisson Sampling
#'
#' @param pik vector of first-order inclusion probabilities
#'
#' @importFrom sampling .as_int UPMEpiktildefrompik UPMEqfromw UPMEsfromq
#'
#'
#' @note this functions is a modified version of function \code{\link[sampling]{UPmaxentropy}},
#' in the \pkg{sampling} package.
#'


pre_CPS <- function(pik){
    #compute matrix q of selection probabilities
    nn = sum(pik)
    nn = sampling::.as_int(nn)
    pik2 = pik[pik != 1]
    nn = sum(pik2)
    nn = sampling::.as_int(nn)
    piktilde = sampling::UPMEpiktildefrompik(pik2)
    w = piktilde/(1 - piktilde)
    q = sampling::UPMEqfromw(w, nn)

    return(q)
}

#' Conditional Poisson Sampling (maximum entropy sampling)
#'
#' Draw a sample by means of Conditional Poisson Sampling
#'
#' @param pik vector of first-order inclusion probabilities
#' @param N population size (excluding self-selecting units)
#' @param q matrix of selection probabilities, computed by means of function
#' \code{UPMEsfromq()} of package \code{sampling}.
#'
#'
#' @note this functions is a modified version of function \code{\link[sampling]{UPmaxentropy}},
#' in the \pkg{sampling} package.
#'
#'

maxEntropy <- function(pik, N, q){
    #select samples
    s2 = sampling::UPMEsfromq(q)
    s = rep(0, times = N)
    s[pik == 1] = 1
    s[pik != 1][s2 == 1] = 1
    return(s)
}
