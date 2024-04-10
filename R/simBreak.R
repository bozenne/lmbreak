### simBreak.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Apr  5 2024 (15:37) 
## Version: 
## Last-Updated: apr 10 2024 (15:52) 
##           By: Brice Ozenne
##     Update #: 41
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * simBreak (documentation)
##' @title Simulate Outcome Follow a Breakpoint Model
##' @description Simulate an outcome variable (Y) related to a covariate (X) according to broken lines
##' @name simBreak
##' 
##' @param n [integer,>0] number of individuals for which a broken line should be simulated.
##' Can contain a second element to indicate the number of observations per individual.
##' @param breakpoint [numeric vector] the position of the breakpoints.
##' Its last element indicates the largest possible X value.
##' @param slope [numeric vector] the regression coefficient relative to each broken line.
##' Should have length the length of argument breakpoint minus one.
##' @param sigma [numeric vector] the residual standard deviation of the outcome relative to each broken line.
##' Should have length the length of argument breakpoint minus one.
##' @param seed [integer, >0] Random number generator (RNG) state used when starting imputation. If NULL no state is set.
##' @param rdist.X [numeric vector or function ] X values for each individuals or function used to generate the X values for each individual.
##' 
##' @keywords datagen
##' 
##' @examples
##' #### different timepoints for each individual ####
##' df1 <- simBreak(10, breakpoint = c(0,1,3,4), slope = c(1,0,-1), sigma = 0.1)
##' if(require(ggplot2)){
##' ggplot(df1, aes(x = X, y = Y, group = id)) + geom_point() + geom_line()
##' }
##'
##' ## with more measurements per individual
##' df1.bis <- simBreak(c(10,50), breakpoint = c(0,1,3,4), slope = c(1,0,-1), sigma = 0.1)
##' if(require(ggplot2)){
##' ggplot(df1.bis, aes(x = X, y = Y, group = id)) + geom_point() + geom_line()
##' }
##' 
##' #### same timepoint for all individuals ###
##' df2 <- simBreak(10, breakpoint = c(0,1,3,4), slope = c(1,0,-1),
##'                 rdist.X = seq(0,4,0.5), sigma = 0.1)
##' if(require(ggplot2)){
##' ggplot(df2, aes(x = X, y = Y, group = id)) + geom_point() + geom_line()
##' }
##' 

## * simBreak (code)
##' @export
simBreak <- function(n, breakpoint, slope, sigma = 1, 
                     rdist.X = stats::runif, seed = NULL){

    ## ** normalize user input
    if(any(n<1)){
        stop("Argument \'n\' should only contain stricly positive integers. \n")
    }
    if(length(n)==1){
        n <- c(n,25)
    }else if(length(n)>2){
        stop("Argument \'n\' should have length at most 2. \n")
    }
    n.breakpoint <- length(breakpoint)
    if(n.breakpoint<1){
        stop("Argument \'breakpoint\' should have length greater or equal to 1. \n")
    }
    if((n.breakpoint-1)!=length(slope)){
        stop("Argument \'slope\' should have the same length as argument \'breakpoint\'. \n")
    }
    
    if((n.breakpoint-1)!=length(sigma)){
        if(length(sigma)==1){
            sigma <- rep(sigma, n.breakpoint-1)
        }else{
            stop("Argument \'sigma\' should have the same length as argument \'breakpoint\'. \n")
        }
    }

    if(!is.character(rdist.X) && !is.function(rdist.X) && !is.numeric(rdist.X)){
        stop("Argument \'rdist.X\' should refer to a function or be a vector of values.")
    }

    ## ** set seed for reproducibility
    if(!is.null(seed)){
        if(!is.null(get0(".Random.seed"))){ ## avoid error when .Random.seed do not exists, e.g. fresh R session with no call to RNG
            old <- .Random.seed # to save the current seed
            on.exit(.Random.seed <<- old) # restore the current seed (before the call to the function)
        }else{
            on.exit(rm(.Random.seed, envir=.GlobalEnv))
        }
        set.seed(seed)        
    }

    ## ** simulate
    if(is.numeric(rdist.X)){
        out <- do.call(rbind,lapply(1:n[1], function(iID){ data.frame(id = iID, X = rdist.X)}))
    }else{
        out <- do.call(rbind,lapply(1:n[1], function(iID){ data.frame(id = iID, X = do.call(rdist.X, args = list(n[2], breakpoint[1], breakpoint[n.breakpoint])))}))
    }
    index <- as.numeric(cut(out$X, breakpoint, include.lowest = TRUE))
    intercept <- cumsum(c(breakpoint[1], diff(breakpoint)*slope))[1:(n.breakpoint-1)]
    out$Y <- intercept[index] + slope[index] * (out$X-breakpoint[index]) + stats::rnorm(NROW(out), mean = 0, sd =  sigma[index])
    
    ## ** export
    return(out)
}


##----------------------------------------------------------------------
### simBreak.R ends here
