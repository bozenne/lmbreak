### fitted.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Apr  8 2024 (09:10) 
## Version: 
## Last-Updated: Apr  8 2024 (09:13) 
##           By: Brice Ozenne
##     Update #: 5
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * fitted.lmbreak
##' @export
fitted.lmbreak <- function(object, ...){

    ## ** check user input
    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }
    
    ## ** export
    return(object$model$fitted)

}

##----------------------------------------------------------------------
### fitted.R ends here
