### plot.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Apr  8 2024 (08:25) 
## Version: 
## Last-Updated: apr 10 2024 (16:41) 
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

## * plot.lmbreak
##' @rdname plot.lmbreak
##' @export
plot.lmbreak <- function(x, ...){
    out <- autoplot.lmbreak(x, ...)
    print(out$plot)
    return(invisible(out))

}

## * plot.mlmbreak 
##' @rdname plot.mlmbreak
##' @export
plot.mlmbreak <- function(x, ...){
    out <- autoplot.mlmbreak(x, ...)
    print(out$plot)
    return(invisible(out))

}

##----------------------------------------------------------------------
### plot.R ends here
