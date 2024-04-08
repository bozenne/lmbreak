### coef.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Apr  6 2024 (12:23) 
## Version: 
## Last-Updated: Apr  8 2024 (11:52) 
##           By: Brice Ozenne
##     Update #: 9
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

##' @export
coef.lmbreak <- function(object, type = "breakpoint", ...){

    ## ** normalize user input
    type <- match.arg(type, c("bp","breakpoint","Us","slope","Vs","discontinuity","lm"))

    ## ** extract coefficient
    if(type %in% c("bp","breakpoint")){
        out <- object$breakpoint$value
    }else if(type %in% c("Us","slope")){
        if(object$continuity == FALSE && !is.null(attr(object$model,"continuity"))){
            out <- coef(attr(object$model,"continuity"))[object$breakpoint$Us]
        }else{
            out <- object$coef[object$breakpoint$Us]
        }
    }else if(type %in% c("Vs","discontinuity")){
        out <- object$coef[object$breakpoint$Vs]
    }else if(type == "lm"){
        otu <- coef(object$model)
    }

    ## ** output
    return(out)
}

##----------------------------------------------------------------------
### coef.R ends here
