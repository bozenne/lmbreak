### plot.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Apr  8 2024 (08:25) 
## Version: 
## Last-Updated: Apr  8 2024 (10:10) 
##           By: Brice Ozenne
##     Update #: 3
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * plot.lmbreak 
##' @export
plot.lmbreak <- function(x, ...){
    out <- autoplot.lmbreak(x, ...)
    print(out$plot)
    return(invisible(out))

}


##----------------------------------------------------------------------
### plot.R ends here
