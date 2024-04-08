### model.frame.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Apr  8 2024 (08:33) 
## Version: 
## Last-Updated: Apr  8 2024 (08:48) 
##           By: Brice Ozenne
##     Update #: 12
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * model.frame.lmbreak
##' @export
model.frame.lmbreak <- function(formula, newdata = NULL, ...){

    ## ** normalize user input
    variable <- formula$breakpoint.var
    if(inherits(formula, "lmbreak")){
        breakpoint <- formula$breakpoint$value
    }else{
        breakpoint <- formula$breakpoint
    }
    n.breakpoint <- length(breakpoint)
    if(is.null(newdata)){
        newdata <- formula$data
    }

    ## ** add missing columns
    if("Us0" %in% names(newdata) == FALSE){
        newdata$Us0 <- newdata[[variable]]
    }

    for(iPoint in 1:n.breakpoint){ ## iPoint <- 1
        newdata[[paste0("Us",iPoint)]] <- pmax(0,newdata[[variable]]-breakpoint[iPoint])
        newdata[[paste0("Vs",iPoint)]] <- -as.numeric(newdata[[variable]]>breakpoint[iPoint])
    }

    ## ** export
    return(newdata)
}

##----------------------------------------------------------------------
### model.frame.R ends here
