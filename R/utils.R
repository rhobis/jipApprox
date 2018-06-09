### Helper functions -----------------------------------------------------------

#' Exclude self-selecting units
#'
#' Exclude self-selecting units and units with probability zero and returns a list
#' with parameters needed to perform sampling
#'
#' @param pik vector of first-order inclusion probabilities
#' @param eps control value for pik
#'
#' @note the code is taken from package \pkg{sampling}
#'
#' @keywords internal


excludeSSU <- function(pik, eps = 1e-06){
    #from sampling package
    if (any(is.na(pik)))
        stop("there are missing values in the pik vector")
    n = sum(pik)
    n = sampling::.as_int(n)
    list = pik > eps & pik < 1 - eps
    pikb <- pik[list]
    N = length(pikb)
    if (N < 1){
        stop("the pik vector has all elements outside of the range [eps,1-eps]")
    } else {
        return( list(pik = pikb,
                     n = sum(pikb),
                     N = length(pikb),
                     s = pik,
                     list = list)
        )
    }
}

#' Save partial results
#'
#' Write joint inclusion probability estimates on a file every \code{by} replications
#'
#' @param iteration integer indicating the iterations the simulation is at
#' @param replications integer indicating the total number of replications
#' for the simulation
#' @param by integer indicating after how many iterations estimates should be
#' written on a file, must be < \code{replications}
#' @param counts matrix with number of occurrences of couple of units up to
#' current replication
#' @param units id of units for which output should be saved
#' @param file_path path of the file to write
#' @param asDataFrame logical, should output be in form of a data frame?
#'
#'
#' @keywords internal

savepartial <- function(iteration,
                        replications,
                        by,
                        counts,
                        units,
                        file_path,
                        asDataFrame){

    #get data ready
    jips <- counts/iteration
    colnames(jips) <- as.character(units)
    rownames(jips) <- colnames(jips)
    if(asDataFrame){
        # Joint inclusion probabilities, from matrix to data frame
        jip_v <- vector()
        for(i in 2:n){
            jip_v <- c(jip_v, jips[(i-1),i:n])
        }
        jips <- data.frame(i=couples[1,],j=couples[2,], pi_ij=jip_v)
    }

    #number of iteration string to append to file name
    if( iteration<1e03){
        niter <- sprintf("%i",iteration)
    } else if(iteration> 999 & iteration<1e04){
        niter <- sprintf("%sK",substr(iteration,1,1))
    }else if(iteration<1e05){
        niter <- sprintf("%sK",substr(iteration,1,2))
    }else if(iteration<1e06){
        niter <- sprintf("%sK",substr(iteration,1,3))
    }else if(iteration<1e07){
        niter <- sprintf("%s.%sM",substr(iteration,1,1), substr(iteration,2,2))
    }else if(iteration<1e08){
        niter <- sprintf("%s.%sM",substr(iteration,1,2), substr(iteration,3,3))
    }else niter <- sprintf("%i",iteration)

    #write output on file
    write.table(jips, file = paste(file_path,'_',niter,'.txt',sep=''))
    Sys.sleep(.1)
}



#' Transform Joint-Inclusion Probability matrix to data.frame

jipMtoDF <- function(jipmat, id=NULL){
    # Joint inclusion probabilities, from matrix to data frame
    # mat = matrix with joint inclusion probabilities
    # id  = optional, id labels
    n <- nrow(jipmat)
    jip_v <- vector()
    for(i in 2:n){
        jip_v <- c(jip_v, jipmat[(i-1),i:n])
    }
    if(!missing(id)){
        couples <- combn(id,2)
    }else if(missing(id) & !is.null(colnames(jipmat))){
        id <- colnames(jipmat)
        couples <- combn(id,2)
    }else couples <- combn(1:n,2)

    return(data.frame(i=couples[1,],j=couples[2,], pi_ij=jip_v))
}


#' Transform Joint-Inclusion Probability data.frame to matrix
#'
#'

#returns an upper triangular matrix
jipDFtoM <- function(pi_ij){
    n <- (1 + sqrt( 1 + 8*length(pi_ij) ) ) / 2
    M <- matrix(0,n,n)
    r <- 1; i <- 1; j <- n-1
    while(r < n){
        M[r,(r+1):n] <- pi_ij[i:j]
        r <- r + 1
        i <- j + 1
        j <- j + (n-r)
    }
    return(M)
}

#returns a symmetric matrix
jipDFtoM2 <- function(pi_ij){
    n <- (1 + sqrt( 1 + 8*length(pi_ij) ) ) / 2
    M <- matrix(0,n,n)
    r <- 1; i <- 1; j <- n-1
    while(r < n){
        M[r,(r+1):n] <- pi_ij[i:j]
        r <- r + 1
        i <- j + 1
        j <- j + (n-r)
    }
    lt <- lower.tri(M)
    M[lt] <- t(M)[lt]
    return(M)
}


