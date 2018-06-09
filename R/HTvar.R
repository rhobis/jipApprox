#' Variance of the Horvitz-Thompson estimator
#'
#' Compute or estimate the variance of the Horvitz-Thompson total estimator
#' by the Horvitz-Thompson or Sen-Yates-Grundy variance estimators.
#'
#' @param y numeric vector representing the variable object of inference
#' @param pikl matrix of second-order (joint) inclusion probabilities; the diagonal
#' must contain the first-order inclusion probabilities.
#' @param sample boolean indicating if sample values are provided.
#' If \code{sample=TRUE}, the function returnss a sample estimate of the variance,
#' while if \code{sample=FALSE}, the Variance is computed over all population units.
#' Default is TRUE.
#' @param method string, indicating if the Horvitz-Thompson (\code{"HT"}) or the
#' Sen-Yates-Grundy (\code{"SYG"}) estimator should be computed.
#'
#'
#' @details
#' The Horvitz-Thompson variance is defined as
#'
#' \deqn{\sum_{i\in U}\sum_{j \in U} \frac{(\pi_{ij} - \pi_i\pi_j)}{\pi_i\pi_j} y_i y_j }{
#' \sum_U \sum_U [ \pi(ij) - \pi(i)\pi(j) ] y(i) y(j) / [ \pi(i)\pi(j) ] }
#'
#' which is estimated by
#'
#' \deqn{\sum_{i\in U}\sum_{j \in U} \frac{(\pi_{ij} - \pi_i\pi_j)}{\pi_i\pi_j\pi_{ij}} y_i y_j }{
#' \sum_s \sum_s [ \pi(ij) - \pi(i)\pi(j) ] y(i) y(j) / [ \pi(i)\pi(j)\pi(ij) ] }
#'
#'
#' The Sen-Yates-Grundy variance is obtained from the Horvitz-Thompson variance
#' by conditioning on the sample size n, and is therefore only appliable to
#' fixed size sampling designs:
#'
#' \deqn{\sum_{i\in U}\sum_{j > i} (\pi_i\pi_j - \pi_{ij}) \Bigl(\frac{y_i}{\pi_i} - \frac{y_j}{\pi_j} \Bigr)^2 }{
#' \sum__U\sum_{j > i} [ \pi(i)\pi(j) - \pi(ij) ] [ y(i)/\pi(i) - y(j)/\pi(j) ]^2 }
#'
#' Its estimator is
#'
#' \deqn{\sum_{i\in U}\sum_{j > i} \frac{(\pi_i\pi_j - \pi_{ij})}{\pi_{ij}} \Bigl(\frac{y_i}{\pi_i} - \frac{y_j}{\pi_j} \Bigr)^2 }{
#' \sum__U\sum_{j > i} [ \pi(i)\pi(j) - \pi(ij) ] [ y(i)/\pi(i) - y(j)/\pi(j) ]^2 / \pi(ij) }
#'
#'
#'
#'
#' @examples
#'
#' ### Generate population data ---
#' N <- 500; n <- 50
#'
#' set.seed(0)
#' x <- rgamma(500, scale=10, shape=5)
#' y <- abs( 2*x + 3.7*sqrt(x) * rnorm(N) )
#'
#' pik  <- n * x/sum(x)
#' pikl <- jip_approx(pik, method='Hajek')
#'
#' ### Dummy sample ---
#' s   <- sample(N, n)
#'
#'
#' ### Compute Variance ---
#' HTvar(y=y, pikl=pikl, sample=FALSE, method="HT")
#' HTvar(y=y, pikl=pikl, sample=FALSE, method="SYG")
#'
#'
#' ### Estimate Variance ---
#' #' HTvar(y=y[s], pikl=pikl[s,s], sample=TRUE, method="HT")
#' #' HTvar(y=y[s], pikl=pikl[s,s], sample=TRUE, method="SYG")
#'
#'
#' @export
#'

HTvar <- function(y, pikl, sample = TRUE, method = "HT") {

    ### Check input ---
    method <- match.arg(method, c("HT", "SYG"))

    if( !is.numeric(pikl) | !is.matrix(pikl) ) stop("Argument 'pikl' must be a square numerical matrix")
    if( !identical( nrow(pikl), ncol(pikl) ) ) stop("Argument 'pikl' is not a square matrix!")
    if( any( is.na(pikl) ) )                   stop("There are NA values in the 'pikl' argument!")
    if( !identical( length(y), nrow(pikl) ) )  stop("Arguments 'y' and 'pikl' have incompatible dimensions!")
    if( !is.numeric(y) )                       stop("Argument 'y' is not numeric!")
    if( any( is.na(y) ) )                      stop("There are NA values in the 'y' argument!")
    if( !is.logical(sample) | is.na(sample) )  stop("Argument 'sample' must be a logical value (TRUE/FALSE)!")


    # message( ifelse(sample, "Estimating ", "Computing "),
    #          ifelse(identical(method, "HT"), "Horvitz-Thompson ", "Sen-Yates-Grundy "),
    #          "variance ...")


    ### Compute/Estimate variance ---
    pik <- diag(pikl)
    pp  <- outer(pik, pik, '*')

    if( identical(method, "HT") ){

        yy <- outer(y, y, '*')
        D  <- (pikl-pp) / pp
        if( sample )
            D <- D / pikl

        v  <- sum( D*yy )

    }else {

        ef  <- y/pik
        yy  <- outer(ef, ef, '-')**2
        D   <- (pp - pikl)
        if( sample ) D <- D/pikl

        v   <- D*yy
        v[!upper.tri(v)] <- 0
        v   <- sum(v)

    }

    ### Return result ---
    return(v)

}

