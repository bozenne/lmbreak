### test-auto-previous-bug.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jul 18 2024 (10:59) 
## Version: 
## Last-Updated: jul 18 2024 (12:05) 
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

if(FALSE){
    library(testthat)
    library(lmbreak)
}

context("Previous bug")

## * from Friederike Holze < friederike.holze@nru.dk > data: Wednesday 17/07/24 15.48.02
## ISSUE: incorrect slope calculation 
## adapted from her example

test_that("Slope calculation", {
    set.seed(10)
    df <- data.frame(score = c(0:10,10:0,rep(0,10)) + rnorm(32,sd = 0.01),
                     value = 0:31)
    e.test <- lmbreak(score ~ 0 + bp(value, pattern = "110"), data = df)
    
    expect_equal(round(model.tables(e.test)$slope,1)[1:3], c(1,-1,0))
    expect_equal(round(model.tables(e.test)$intercept,1)[1:3], c(0,10.5,0))
})
##----------------------------------------------------------------------
### test-auto-previous-bug.R ends here
