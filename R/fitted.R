### fitted.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Apr  8 2024 (09:10) 
## Version: 
## Last-Updated: jul 18 2024 (10:32) 
##           By: Brice Ozenne
##     Update #: 6
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
fitted.lmbreak <- function(object, continuity = NULL, ...){

    ## ** check user input
    ## *** dots
    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }
    
    ## *** continuity
    if(is.null(continuity)){
        continuity <- (object$opt$continuity==FALSE)
    }

    ## ** extract from object
    if(continuity && !is.null(attr(object$model,"continuity"))){
        model <- attr(object$model,"continuity")
    }else{
        model <- object$model
    }

    out <- model$fitted
    
    ## ** export
    return(out)

}

##----------------------------------------------------------------------
### fitted.R ends here
