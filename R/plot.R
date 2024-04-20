### plot.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Apr  8 2024 (08:25) 
## Version: 
## Last-Updated: Apr 20 2024 (12:17) 
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

## * plot.lmbreak
##' @rdname plot.lmbreak
##' @export
plot.lmbreak <- function(x, y = NULL, ...){
    out <- autoplot.lmbreak(x, y = y, ...)
    print(out$plot)
    return(invisible(out))

}

## * plot.mlmbreak 
##' @rdname plot.mlmbreak
##' @export
plot.mlmbreak <- function(x, y = NULL, ...){
    out <- autoplot.mlmbreak(x, y = y, ...)
    print(out$plot)
    return(invisible(out))

}

##----------------------------------------------------------------------
### plot.R ends here
