### Helper functions -----------------------------------------------------------



#' Transform a matrix of Joint-Inclusion Probabilities to a data.frame
#'
#' @param jip a square matrix of joint-inclusion probabilities, symmetric or
#' upper-triangular
#' @param id optional, vector of id labels, its length should be equal to
#' \code{ncol(jip)} and \code{nrow(jip)}
#'
#'
#' @export

jipMtoDF <- function(jip, id=NULL){

    n <- nrow(jip)
    jip_v <- vector()
    for(i in 2:n){
        jip_v <- c(jip_v, jip[(i-1),i:n])
    }
    if(!missing(id)){
        couples <- combn(id,2)
    }else if(missing(id) & !is.null(colnames(jip))){
        id <- colnames(jip)
        couples <- combn(id,2)
    }else couples <- combn(1:n,2)

    return(data.frame(pi=couples[1,],pj=couples[2,], pij=jip_v))
}


#' Transform a Joint-Inclusion Probability data.frame to a matrix
#'
#' @param jip vector or data.frame containing the joint-inclusion probabilities
#' @param symmetric boolean, if \code{TRUE}, returns a symmetric matrix, otherwise,
#' an upper triangular matrix
#'
#' @return a symmetric matrix of joint-inclusion probabilities if \code{TRUE}, otherwise,
#' an upper triangular matrix
#'
#' @export


jipDFtoM <- function(jip, symmetric = TRUE){
    n <- (1 + sqrt( 1 + 8*length(jip) ) ) / 2
    M <- matrix(0,n,n)
    r <- 1; i <- 1; j <- n-1
    while(r < n){
        M[r,(r+1):n] <- jip[i:j]
        r <- r + 1
        i <- j + 1
        j <- j + (n-r)
    }
    if( symmetric ){
        lt <- lower.tri(M)
        M[lt] <- t(M)[lt]
    }
    return(M)
}




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
#' @param design_name string indicating the name of the sampling design to include
#' in the filename
#' @param counts matrix with number of occurrences of couple of units up to
#' current replication
#' @param units id of units for which output should be saved
#' @param filename name of the file to write on disk
#' @param path path to the directory where the file should be saved
#' @param status \code{1} if partial results are written before the maximum number
#' of replications is reached, \code{0} otherwise
#' @param as_data_frame logical, should output be in form of a data frame?
#'
#'
#' @keywords internal

save_output <- function(iteration,
                        design_name,
                        counts,
                        units,
                        filename,
                        path,
                        status,
                        as_data_frame){

    if( !(status %in% c(0,1)) ) stop("status should be either 0 or 1")
    #get data ready
    M <- ifelse( identical(design_name, "systematic"), iteration, iteration+1)
    jips <- counts/M
    colnames(jips) <- rownames(jips)  <- as.character(units)
    if(as_data_frame) jips <- jipMtoDF(jips, id=units)

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

    #define filename
    if( is.null(filename) ){
        filename <- paste0( file.path(path, design_name), '_', niter,'.txt')
    }else{
        if( as.logical(status) ){
            filename <- paste0(path, '/', niter, '_', filename)
        }else filename <- paste0( file.path(path, filename) )
    }

    #write output on file
    write.table(jips, file = filename)

}



#' Check if a number is integer
#'
#' Check if \code{x} is an integer number, differently from \code{is.integer},
#' which checks the type of the object \code{x}
#'
#' @param x a scalar or a numeric vector
#' @param tol a scalar, indicating the tolerance
#'
#'
#' @note From the help page of function \code{\link[base]{is.integer}}
#'
#' @keywords internal

is.wholenumber <- function(x, tol = .Machine$double.eps^0.5){
    abs(x - round(x)) < tol
}
