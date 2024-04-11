### 0-onload.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr 11 2024 (09:10) 
## Version: 
## Last-Updated: apr 11 2024 (09:19) 
##           By: Brice Ozenne
##     Update #: 4
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

lmbreak.env <- new.env() # create a specific environment for the package

.onAttach <- function(lib, pkg="lmbreak") {
    desc <- utils::packageDescription(pkg)
    packageStartupMessage(desc$Package, " version ",desc$Version)
    lmbreak.options(reinitialise = TRUE) # generate .lmbreak-options when loading the package   
}


##----------------------------------------------------------------------
### 0-onload.R ends here
