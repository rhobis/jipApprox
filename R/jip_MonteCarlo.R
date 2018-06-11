#' Approximate inclusion probabilities by Monte Carlo simulation
#'
#' Approximate first and second-order inclusion probabilities by means of
#' Monte Carlo simulation.
#' Estimates are obtained as proportion of the number of occurrences of each unit or
#' couple of units over the total number of replications. One unit is added to both
#' numerator and  denominator to assure strict positivity of estimates (Fattorini, 2006).
#'
#' @param x size measure or first-order inclusion probabilities
#' @param n sample size (for fixed size samplings), or expected sample size (for Poisson sampling)
#' @param replications number of independent Monte Carlo replications
#' @param design sampling procedure to be used for sample selection.
#'        Either a string indicating the name of the sampling design or a function;
#'        see section "Details" for more information.
#' @param units id of units for which probabilities have to be estimated.
#'        Optional, if missing, estimates are produced for the whole population
#' @param seed a valid seed value for reproducibility
#' @param as_data_frame logical, should output be in a data.frame form? if FALSE,
#' a matrix is returned
#' @param design_pars only used when a function is passed to argument \code{design},
#'        list of parameters to pass to the sampling design function.
#' @param write_on_file logical, should output be written on a text file?
#' @param file_path string indicating the path to the directory where the text
#' file with the output should be created, used only if \code{write_on_file} is TRUE.
#' It should never end with a trailing slash '/',
#' @param by optional; integer scalar indicating every how many replications a partial output
#' should be saved
#'
#' @details
#' Argument \code{design} accepts either a string indicating the sampling design
#' to use to draw samples or a function.
#' Accepted designs are "brewer", "tille", "maxEntropy", "poisson",
#' "sampford", "systematic", "randomSystematic", "srs", chao", "sunter", "srs", "fpdust", "fpdust_pps".
#' The user may also pass a function as argument; such function should take as input
#' the parameters passed to argument \code{design_pars} and return either a logical
#' vector or a vector of 0 and 1,  where \code{TRUE} or \code{1} indicate sampled
#' units and \code{FALSE} or \code{0} indicate non-sample units.
#' The length of such vector must be equal to the length of \code{x}
#' if \code{units} is not specified, otherwise it must have the same length of \code{units}.
#'
#'
#' @return A matrix of estimated inclusion probabilities if \code{as_data_frame=FALSE},
#' otherwise a data.frame with three columns: the first two indicate the id of
#' the couple of units, while the third one contains the joint-inclusion probability
#' values. Please, note that when \code{as_data_frame=TRUE}, first-order
#' inclusion probabilities are not returned.
#'
#'
#' @examples
#' ### Generate population data ---
#' N <- 20; n<-5
#'
#' set.seed(0)
#' x <- rgamma(N, scale=10, shape=5)
#' y <- abs( 2*x + 3.7*sqrt(x) * rnorm(N) )
#'
#' pik  <- n * x/sum(x)
#'
#' ### Approximate joint-inclusion probabilities
#' pikl <- jip_MonteCarlo(x=pik, n = n, replications = 100, design = "brewer")
#' pikl <- jip_MonteCarlo(x=pik, n = n, replications = 100, design = "tille")
#' pikl <- jip_MonteCarlo(x=pik, n = n, replications = 100, design = "maxEntropy")
#' pikl <- jip_MonteCarlo(x=pik, n = n, replications = 100, design = "randomSystematic")
#' pikl <- jip_MonteCarlo(x=pik, n = n, replications = 100, design = "systematic")
#' pikl <- jip_MonteCarlo(x=pik, n = n, replications = 100, design = "sampford")
#' pikl <- jip_MonteCarlo(x=pik, n = n, replications = 100, design = "poisson")
#'
#'
#'
#' @export
#'
#' @importFrom sampling inclusionprobabilities srswor UPrandomsystematic UPpoisson
#'

jip_MonteCarlo <- function(x, n, replications = 1e06,
                           design,
                           units,
                           seed = NULL,
                           as_data_frame = FALSE,
                           design_pars,
                           write_on_file = FALSE,
                           file_path,
                           by = NULL
){

    ### Check input ------------------------------------------------------------

    design <- match.arg(design, c('brewer',
                                  'tille',
                                  'maxEntropy',
                                  'randomSystematic',
                                  'sampford',
                                  'poisson',
                                  'systematic')
    )

    ## x
    if (is.data.frame(x))
        if (ncol(x) > 1)
            stop("Argument x is not a vector")
    else X = unlist(x)
    else if (is.matrix(x))
        if (ncol(x) > 1)
            stop("Argument x is not a vector")
    else X = X[, 1]
    else if (is.list(x))
        if (length(x) > 1)
            stop("Argument x is not a vector")
    else X = unlist(x)

    if (any(is.na(x))){
        stop("Argument x has some missing values")
    }
    ## units
    if(missing(units)){
        units <- seq_along(x)
    }else if( is.null(units) ){
        units <- seq_along(x)
    }else if(length(units) < 2){
        stop("Not a valid units vector, it should include at least two unit ids")
    } else if(!is.numeric(units)){
        stop("units is not a numeric vector")
    }
    ##n
    if(!is.numeric(n)){
        stop("n must be a numeric value")
    }else n <- sampling::.as_int(ceiling(n))


    # seed, check if it's numeric
    if( !missing(seed) && !is.numeric(seed)){
        stop("Please, provide a numeric seed")
    }
    ## replications
    if(!is.numeric(replications) || replications < 2){
        stop("Not a valid replications number, it should be a number greater than 1")
    } else if(replications < 1000)
        message("The number of Monte Carlo replications is too small, estimates may be unreliable!")



    ## file_path
    if( write_on_file ){
        #file_path
        if( missing(file_path) ) file_path <- getwd()
        if( !is.character(file_path) ) stop("file_path should be a string!")
        if(!dir.exists(file_path)){
            dc <- dir.create(file_path, showWarnings = TRUE, recursive = TRUE)
            if( !dc ) stop("There is no folder at the path specified in file_path ",
                           "and the creation of a new directory failed, please provide a path to an existing directory!")
        }
        #by
        if( is.null(by) ) by <- replications
        if( !is.numeric(by) ) stop("Argument 'by' should be a positive integer scalar")
        if( by > replications) by <- replications

        # print path
        message("Output will be written in the directory: ", file_path)

    }



    ### Initialisation ----
    pik <- sampling::inclusionprobabilities(x,n)
    k <- length(units)
    if( identical(design, "systematic") ){
        counts <- matrix(0,k,k)
    }else counts <- matrix(1,k,k)

    pb <- txtProgressBar(min = 0, max = replications, style = 3) # Progress bar

    if( is.character(design)){
        ### Preliminary computations (for some sampling methods) ---
        pre <- switch(EXPR=design,
                      'default' = NULL,
                      'brewer' = excludeSSU(pik),
                      'sampford' = excludeSSU(pik),
                      'tille' = pre_tille(pik=pik),
                      'maxEntropy' = pre_CPS(pik=pik)
        )

        #parameters for sampling function
        pars <- switch(EXPR=design,
                       'brewer' = pre,
                       'tille' = pre,
                       'maxEntropy' = list(pik=pik, N=length(pik), q=pre),
                       'randomSystematic' = list(pik=pik),
                       'sampford' = pre,
                       'poisson' = list(pik=pik),
                       'systematic' = list(pik=pik)
        )
        # sampling function
        smplFUN <- switch(EXPR=design,
                          'brewer' = brewer,
                          'tille' = tille,
                          'maxEntropy' = maxEntropy,
                          'randomSystematic' = sampling::UPrandomsystematic,
                          'sampford' = sampford,
                          'poisson' = sampling::UPpoisson,
                          'systematic' = sampling::UPsystematic
        )
    }else if( is.function(design) ){
        pars    <- design_pars
        smplFUN <- design
    }else stop("Argument design is not well-specified: it should be either a string representing ",
               "one of the available sampling designs or an object of class function!")

    ### Monte Carlo simulation ----
    set.seed(seed)
    for(r in 1:replications){
        setTxtProgressBar(pb, r)

        s <- do.call(smplFUN,pars) ### vector of 0 and 1
        counts <- counts + outer(s[units],s[units], "*")
        if(write_on_file){
            if(r %% by == 0){
                savepartial(r, design, counts, units, file_path, as_data_frame)
            }
        }
    }
    close(pb)
    ### ----

    ### Return estimates ----
    M <- ifelse( identical(design, "systematic"), replications, replications+1)
    jips <- counts/M
    colnames(jips) <- rownames(jips) <- as.character(units)

    # change class, from matrix to data frame
    if(as_data_frame) jips <- jipMtoDF(jips, id=units)

    ### Output
    return(jips)

}
