### lmbreak.options.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Apr 16 2021 (12:01) 
## Version: 
## Last-Updated: apr 12 2024 (17:17) 
##           By: Brice Ozenne
##     Update #: 177
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * lmbreak.options (documentation) 
#' @title Global options for lmbreak package
#' @include 0-onload.R
#'
#' @description Update or select global options for the lmbreak package.
#'
#' @param ... options to be selected or updated
#' @param reinitialise should all the global parameters be set to their default value
#'
#' @details The options are: \itemize{
#' \item enforce.continuity [logical]: enforce.continuity [logical] in the case where no continuous solution could be found,
#' should a non-continuous breakpoint model be kept (\code{FALSE}) or a continuous breakpoint model be enforced (\code{TRUE}, refiting without the Vs terms). Used by \code{lmbreak}.
#' \item optimizer [character]: breakpoints are estimated using direct minimisation of the residual sum of squares (\code{"optim"}) or using the approach proposed by Muggeo et al. in 2003 (\code{"Muggeo"}).
#' \item init.gam [logical]: should a spline model be used to define initialization points (where the first derivative of the spline changes sign).
#' \item init.quantile [integer]: number of quantiles used to try different initializations.
#' \item minR2 [numeric,0-1]: minimum R2 when selecting the pattern for qualifying the convergence as satisfactory.
#' \item n.iter [integer, >0] the maximum number of iterations used to estimates the breakpoints. Used by \code{lmbreak}.
#' \item optimize.step [logical]: should a full update of the breakpoint value be performed at each step (\code{FALSE})
#' or a partial update minimizing the sum of residual absolute value of the model without Vs terms.
#' \item tol [numeric, >0] the maximum accpetable difference between two consecutive estimates of the breakpoints. Used by \code{lmbreak}.
#' Can have an optional second element used to assess the continuity of the estimated fit at the breakpoint (i.e. maximal acceptable magnitude of the Vs terms).
#' }
#'
#' @return A list containing the default options.
#'
#' @references Muggeo, V. M. R. Estimating regression models with unknown break-points.
#' Statistics in medicine 2003; 22:3055-3071.
#' 
#' @keywords utilities


 
## * lmbreak.options (code)
#' @export
lmbreak.options <- function(..., reinitialise = FALSE){

    default <- list(enforce.continuity = TRUE,
                    init.gam = TRUE,
                    init.quantile = 0.5,
                    optimizer = "Muggeo",
                    minR2 = 0.01,
                    n.iter = 20,
                    optimize.step = 0,
                    tol = c(1e-3,1e-2))

    if (reinitialise == TRUE) {
        assign(".lmbreak-options", value = default, envir = lmbreak.env)
    
    return(invisible(get(".lmbreak-options", envir = lmbreak.env)))
    
  }else{
    
      args <- list(...)

      if(!is.null(names(args))){
          object <- get(".lmbreak-options", envir = lmbreak.env)
      }else{
          object <- try(get(".lmbreak-options", envir = lmbreak.env))
          if(inherits(object,"try-error")){
              object <- default
          }
      }

      if(length(args)==0){ ## read all
          return(object)
      }else if (!is.null(names(args))) { ## write

          if(any(names(args) %in% names(object) == FALSE)){
              stop("Incorrect element selected: \"",paste0(names(args)[names(args) %in% names(object) == FALSE], collapse = "\" \""),"\"\n",
                   "Available elements: \"",paste0(setdiff(names(object),names(args)), collapse = "\" \""),"\"\n")
          }

          if("n.iter" %in% names(args) && (!is.numeric(args$n.iter) || args$n.iter <= 0 || args$n.iter %% 1 > 0) ){
              stop("Argument \'n.iter\' must be a strictly positive integer. \n")
          }
          if("tol" %in% names(args) && (!is.numeric(args$tol) || args$tol <= 0)){
              stop("Argument \'tol\' must be a strictly positive numeric value. \n")
          }
          if("tol" %in% names(args) && length(args$tol)!=2){
              stop("Argument \'tol\' must be have length 2. \n")
          }
          if("enforce.continuity" %in% names(args) && !is.logical(args$enforce.continuity)){
              stop("Argument \'enforce.continuity\' must be of type logical. \n")
          }
          if("optimize.step" %in% names(args) && !is.logical(args$optimize.step)){
              stop("Argument \'optimize.step\' must be of type logical. \n")
          }
          if("init.gam" %in% names(args) && !is.logical(args$init.gam)){
              stop("Argument \'init.gam\' must be of type logical. \n")
          }
          object[names(args)] <- args
      
          assign(".lmbreak-options", 
                 object, 
                 envir = lmbreak.env)
      
          return(invisible(object))
      
      } else {# read
          args <- unlist(args)
          if(any(args %in% names(object) == FALSE)){
              stop("Incorrect element selected: \"",paste0(args[args %in% names(object) == FALSE], collapse = "\" \""),"\"\n",
                   "Available elements: \"",paste0(setdiff(names(object),args), collapse = "\" \""),"\"\n")
          }
          return(object[args])
      }
    
  }
  
  
}


##----------------------------------------------------------------------
### lmbreak.options.R ends here
