### lmbreak.options.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Apr 16 2021 (12:01) 
## Version: 
## Last-Updated: apr 11 2024 (09:20) 
##           By: Brice Ozenne
##     Update #: 159
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
#' \item n.iter [integer, >0] the maximum number of iterations used to estimates the breakpoints. Used by \code{lmbreak}.
#' \item tol [numeric, >0] the maximum accpetable difference between two consecutive estimates of the breakpoints. Used by \code{lmbreak}.
#' \item enforce.continuity [logical]: enforce.continuity [logical] in the case where no continuous solution could be found,
#' should a non-continuous breakpoint model be kept (\code{FALSE}) or a continuous breakpoint model be enforced (\code{TRUE}, refiting without the Vs terms). Used by \code{lmbreak}.
#' }
#'
#' @return A list containing the default options.
#'
#' @keywords utilities


 
## * lmbreak.options (code)
#' @export
lmbreak.options <- function(..., reinitialise = FALSE){

    default <- list(n.iter = 50,
                    tol = 1e-3,
                    enforce.continuity = TRUE)

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
              stop("Argument \'enforce.continuity\' must be a strictly positive numeric value. \n")
          }
          if("enforce.continuity" %in% names(args) && !is.logical(args$enforce.continuity)){
              stop("Argument \'enforce.continuity\' must be of type logical. \n")
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
