#' Approximate Joint-Inclusion Probabilities
#'
#' Approximations of joint-inclusion probabilities by means of first-order
#' inclusion probabilities.
#'
#' @param pik numeric, vector of first-order inclusion probabilities for all
#' population units.
#' @param method string representing one of the available approximation methods.
#'
#'
#' @details
#' Available methods are \code{"Hajek"}, \code{"HartleyRao"}, \code{"Tille"},
#' \code{"Brewer1"},\code{"Brewer2"},\code{"Brewer3"}, and \code{"Brewer4"}.
#' Note that these methods were derived for high-entropy sampling designs.
#'
#' Hájek (1964) approximation is derived under Maximum Entropy sampling design
#' and is given by
#'
#' \deqn{\tilde{\pi}_{ij} = \pi_i\pi_j \frac{\bigl( 1 - (1-\pi_i)(1-\pi_j)}{d} \bigr)}{
#'       \pi(ij) = \pi(i) \pi(j) (1 - ( 1-\pi(i) )( 1 -\pi(j) ) )/d}
#'  where \eqn{d = \sum_{i\in U} \pi_i(1-\pi_i) }{d = \sum \pi(i)(1-\pi(i))}
#'
#' Hartley and Rao (1962) proposed the following approximation under
#' randomised systematic sampling:
#'
#' \deqn{\tilde{\pi}_{ij} = \frac{n-1}{n} \pi_i\pi_j + \frac{n-1}{n^2} (\pi_i^2 \pi_j + \pi_i \pi_j^2)
#'       - \frac{n-1}{n^3}\pi_i\pi_j \sum_{i\in U} \pi_j^2 + \frac{2(n-1)}{n^3} (\pi_i^3 \pi_j + \pi_i\pi_j^3 + \pi_i^2 \pi_j^2)
#'       - \frac{3(n-1)}{n^4} (\pi_i^2 \pi_j + \pi_i\pi_j^2) \sum_{i \in U}\pi_i^2
#'       + \frac{3(n-1)}{n^5} \pi_i\pi_j \biggl( \sum_{i\in U} \pi_i^2 \biggr)^2
#'       - frac{2(n-1)}{n^4} \pi_i\pi_j \sum{i \in U} \pi_j^3  }{*see pdf version of documentation*}
#'
#' Tillé (1996) proposed the approximation \eqn{\tilde{\pi}_{ij} = \beta_i\beta_j}{\pi(ij) = \beta_i \beta_j},
#' where the coefficients \eqn{\beta_i}{\beta} are computed iteratively through the
#'    following procedure:
#'     \enumerate{
#'         \item \eqn{\beta_i^{(0)} = \pi_i, \,\, \forall i\in U}{\beta(0) = \pi, i = 1, ..., N}
#'         \item \eqn{ \beta_i^{(2k-1)} = \frac{(n-1)\pi_i}{\beta^{(2k-2)} - \beta_i^{(2k-2)}}  }{
#'                     \beta(2k-1) = ( (n-1)\pi )/(\sum\beta(2k-2) - \beta(2k-2)) }
#'         \item \eqn{\beta_i^{2k} = \beta_i^{(2k-1)}
#'         \Biggl( \frac{n(n-1)}{(\beta^(2k-1))^2 - \sum_{i\in U} (\beta_k^{(2k-1)})^2 } \Biggr)^(1/2) }{
#'         \beta(2k) = \beta(2k-1) ( n(n-1) / ( (\sum\beta(2k-1))^2 - \sum( \beta(2k-1)^2 ) ) )^(1/2) }
#'     }
#'     \eqn{  \text{with} \beta^{(k)} = \sum_{i\in U} \beta_i^{i}, \,\, k=1,2,3, \dots }{}
#'
#' Finally, Brewer (2002) and Brewer and Donadio (2003) proposed four approximations,
#' which are defined by the general form
#'
#' \deqn{\tilde{\pi}_{ij} = \pi_i\pi_j (c_i + c_j)/2  }{ \pi(ij) = \pi(i)\pi(j) [c(i) + c(j) ]/2 }
#'
#' where the \eqn{c_i} determine the approximation used:
#'
#' \itemize{
#'     \item Equation (9), \code{method="Brewer1"}:
#'         \deqn{c_i = (n-1) / (n-\pi_i)}{c(i) = [n-1] / [n-\pi(i) ]}
#'    \item Equation (10), \code{method="Brewer2"}:
#'         \deqn{c_i = (n-1) / (n- n^{-1}\sum_{i\in U}\pi_i^2)}{c(i) = [n-1] / [n- \sum_U \pi(i)^2 / n ]}
#'     \item Equation (11), \code{method="Brewer3"}:
#'         \deqn{c_i = (n-1) / (n - 2\pi_i + n^{-1}\sum_{i\in U}\pi_i^2)}{c(i) = [n-1] / [n- 2\pi(i) + \sum_U \pi(i)^2 / n ]}
#'     \item Equation (18), \code{method="Brewer4"}:
#'         \deqn{c_i = (n-1) / (n - (2n-1)(n-1)^{-1}\pi_i + (n-1)^{-1}\sum_{i\in U}\pi_i^2)}{
#'         c(i) = [n-1] / [n- \pi(i)(2n -1)/(n-1) + \sum_U \pi(i)^2 / (n-1) ]}
#'
#' }
#'
#'
#'
#'
#' @return A symmetric matrix of inclusion probabilities, which diagonal is the
#' vector of first-order inclusion probabilities.
#'
#'
#' @references
#'
#' Hartley, H.O.; Rao, J.N.K., 1962. Sampling With Unequal Probability and Without Replacement.
#' The Annals of Mathematical Statistics 33 (2), 350-374.
#'
#' Hájek, J., 1964. Asymptotic Theory of Rejective Sampling with Varying Probabilities from a Finite Population.
#' The Annals of Mathematical Statistics 35 (4), 1491-1523.
#'
#' Tillé, Y., 1996. Some Remarks on Unequal Probability Sampling Designs Without Replacement.
#' Annals of Economics and Statistics 44, 177-189.
#'
#' Brewer, K.R.W.; Donadio, M.E., 2003. The High Entropy Variance of the Horvitz-Thompson Estimator.
#' Survey Methodology 29 (2), 189-196.
#'
#'
#'
#' @examples
#'
#'### Generate population data ---
#' N <- 20; n<-5
#'
#' set.seed(0)
#' x <- rgamma(N, scale=10, shape=5)
#' y <- abs( 2*x + 3.7*sqrt(x) * rnorm(N) )
#'
#' pik  <- n * x/sum(x)
#'
#' ### Approximate joint-inclusion probabilities ---
#' pikl <- jip_approx(pik, method='Hajek')
#' pikl <- jip_approx(pik, method='HartleyRao')
#' pikl <- jip_approx(pik, method='Tille')
#' pikl <- jip_approx(pik, method='Brewer1')
#' pikl <- jip_approx(pik, method='Brewer2')
#' pikl <- jip_approx(pik, method='Brewer3')
#' pikl <- jip_approx(pik, method='Brewer4')
#'
#'
#'
#' @export



jip_approx <- function( pik, method ){

    ### Check input ---
    method <- match.arg(method, c( "Hajek",
                                   "HartleyRao",
                                   "Tille",
                                   "Brewer1",
                                   "Brewer2",
                                   "Brewer3",
                                   "Brewer4")
    )

    if( !identical( class(pik), "numeric" ) ){
        stop( "Argument 'pik' should be a numeric vector!")
    }else if( length(pik) < 2 ){
        stop( "The 'pik' vector is too short!" )
    }else if( any(pik<0)  | any(pik>1) ){
        stop( "Some values of the 'pik' vector are outside the interval [0, 1]")
    }else if( !isTRUE(all.equal(sum(pik), as.integer(sum(pik)) )) ){
        stop( "The sum of 'pik' values is not an integer!")
    }

    ### Call method ---
    jips <- switch(method,
                   "Hajek"       = jip_Hajek(pik),
                   "HartleyRao"  = jip_HartleyRao(pik),
                   "Tille"       = jip_Tille(pik),
                   "Brewer1"    = jip_Brewer(pik, method),
                   "Brewer2"    = jip_Brewer(pik, method),
                   "Brewer3"    = jip_Brewer(pik, method),
                   "Brewer4"    = jip_Brewer(pik, method)
    )

    ### Return result ---
    return( jips )
}



#' Brewer's joint-inclusion probability approximations
#'
#' Approximation of joint inclusion probabilities by one of the estimators
#' proposed by Brewer and Donadio (2003)
#'
#' @inheritParams jip_approx
#'
#' @details
#' \code{"Brewer18"} is the approximation showed in equation (18) of Brewer and Donadio (2003)
#'

jip_Brewer <- function(pik, method){

    method <- match.arg(method, c("Brewer1",
                                  "Brewer2",
                                  "Brewer3",
                                  "Brewer4")
    )
    n <- sum(pik)

    ### Compute c values ---
    if( identical(method, "Brewer1") ){ #Equation (9)

        ci <- (n-1) / (n-pik)

    }else if( identical(method, "Brewer2") ){ #Equation (10)

        ci <- (n-1) / (n - sum(pik**2)/n)

    }else if( identical(method, "Brewer3") ){ #Equation (11)

        ci <- (n-1) / (n - 2*pik + sum(pik**2)/n)

    }else if( identical(method, "Brewer4") ){ #Equation (18)

        ci <- (n-1) / ( n - (2*n-1)/(n-1)*pik + sum(pik**2)/(n-1) )

    }

    ### Estimate jips ---
    if( identical(method, "Brewer2") ){
        cc <- ci*2
    }else cc <- outer(ci,ci, '+')
    pp <- outer(pik,pik, '*')
    out <- pp*cc / 2
    diag(out) <- pik

    ### Return result ---
    return(out)
}


#' Hájek's joint-inclusion probability approximation
#'
#' Estimate joint-inclusion probabilities using Hájek (1964) equation
#'
#'
#' @inheritParams jip_approx
#'
#'

jip_Hajek <- function(pik){

    ### Estimate jips ---
    d   <- sum(pik*(1-pik))
    out <- outer(pik,pik,'*') * (1 - outer(1-pik, 1-pik, '*') / d)
    diag(out) <- pik

    ### Return result ---
    return(out)
}




#' Hartley-Rao approximation of joint-inclusion probabilities
#'
#' Approximation of joint-inclusion probabilities with precision of order
#' \eqn{O(N^{-4})} for the random systematic sampling design
#' by Hartley and Rao (1962), pag. 369 eq. 5.15
#'
#' @inheritParams jip_approx
#'

jip_HartleyRao <- function(pik){

    ### Estimate jips ---
    pp <- outer(pik, pik, '*')
    n  <- sum(pik)
    nn <- (n-1)/n

    p2 <- pik**2
    p3 <- pik**3
    sp2 <- sum(p2)
    sp3 <- sum(p3)

    out <-  nn*pp +
        nn/n * (outer(p2,pik) + outer(pik, p2)) -
        nn/(n**2) * pp * sp2 +
        (2*nn)/(n**2) * (outer(p3,pik) + outer(pik,p3) + outer(p2,p2)) -
        3*nn/n**3 * (outer(p2,pik) + outer(pik,p2)) * sp2 +
        3*nn/n**4 * pp * sp2**2 -
        2*nn/n**3 * pp * sp3

    diag(out) <- pik

    ### Return result ---
    return(out)
}



#' Tillé's approximation of joint-inclusion probabilities
#'
#' Compute the approximation of joint-inclusion probabilities by means of the
#' Iterative Proportional Fitting Procedure (IPFP) proposed by Tillé (1996)
#'
#' @param maxIter a scalar indicating the maximum number of iterations for the
#' fixed-point procedure
#' @param eps tolerance value for the convergence of the fixed-point procedure
#' @inheritParams jip_approx
#'

jip_Tille <- function(pik,eps=1e-06, maxIter=1000){

    n <- sum(pik)

    ### Compute b_k values ---
    b0   <- pik
    iter <- 0
    err  <- Inf
    while( iter<maxIter & err>eps ){
        b1 <- (n-1)*pik / ( sum(b0)-b0 )
        b2 <- b1 * sqrt( n*(n-1) / ( sum(b1)**2 - sum(b1**2) ) )
        b0 <- b2
        tab <- outer(b2,b2,'*');    diag(tab) <- 0
        err <- max( abs( rowSums(tab) - pik*(n-1) ) )
        iter <- iter + 1
    }

    ### Estimate jips ---
    out <- outer(b2,b2,'*')
    diag(out) <- pik


    ### Return result ---
    return(out)
}


