#' @title Tillé's elimination procedure - elimination probabilities
#'
#' @description Computes a matrix with elimination probabilities for each step
#' of Tillé's elimination procedure
#'
#'
#' @param pik vector of first-order inclusion probabilities

pre_tille <- function(pik){ #computations needed just one time
    out <- excludeSSU(pik)

    n <- out$n
    N <- out$N
    pik <- out$pik

    b = rep(1, N)
    #matrix of elimination probabilities
    mv <- matrix(0, nrow=N-n, ncol=N)
    for (i in 1:(N - n)){
        a = sampling::inclusionprobabilities(pik, N - i)
        v = 1 - a/b
        b = a
        mv[i,] = v
    }

    out$pmat <- mv
    return(out)
}


#' @title Tillé's elimination procedure
#'
#' @description Draw a sample by means of Tillé's elimination procedure
#'
#'
#' @param pik vector of first-order inclusion probabilities
#' @param n sample size (excluding self-selecting units)
#' @param N population size (excluding self-selecting units)
#' @param s vector of length N, with 1s at the positions of self-selecting units
#' @param list vector with positions of self selcting units
#' @param pmat matrix of dimension $(N-n)xN, where each row has elimination probabilities
#' for population units for one step of the procedure.
#'
#'
#' @export

tille <- function(pik, n, N, s, list, pmat){ #sample selection
    sb = rep(1, N)
    for(j in 1:(N-n)){
        p = pmat[j,] * sb
        p = cumsum(p)
        u = runif(1)
        lp <- length(p)
        for (h in 1:lp){
            if (u < p[h]) break
        }
        sb[h] = 0
    }
    s[list] <- sb
    return(s)
}
