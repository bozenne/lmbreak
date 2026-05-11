### test-auto-previous-bug.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jul 18 2024 (10:59) 
## Version: 
## Last-Updated: maj 11 2026 (17:16) 
##           By: Brice Ozenne
##     Update #: 8
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

## * from nowahah 21/04/26
## handling error in optim
## Error in stats::optim(par = initialization.trans, method = "BFGS", fn = function(psi) { : 
##   non-finite value supplied by optim
test_that("Handling optim error", {
    data.test <- data.frame("ID" = c("5", "5", "5", "5", "5", "5", "5", "5", "6", "6", "6", "6", "6", "6", "6", "6", "7", "7", "7", "7", "7", "7", "7", "7", "7", "7"), 
                            "time" = c(  0,  23,  41,  60,  79,  99, 118, 141,   0,  16,  37,  58,  78, 101, 119, 140,   0,  20,  40,  58,  77, 101, 122, 139, 159, 178), 
                            "segment" = c("1", "1", "1", "1", "2", "2", "2", "2", "1", "1", "1", "1", "2", "2", "2", "2", "1", "1", "1", "1", "1", "1", "2", "2", "2", "2"), 
                            "pattern" = c(1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0), 
                            "true.traj" = c(0.000000, 3.049641, 5.436317, 7.955586, 9.133016, 9.133016, 9.133016, 9.133016, 0.000000, 2.051788, 4.744760, 7.437731, 9.654467, 9.654467, 9.654467, 9.654467, 0.000000, 1.685005, 3.370010, 4.886514, 6.487269, 8.509275, 8.670334, 8.670334, 8.670334, 8.670334), 
                            "noised.traj" = c(-0.1739980,  2.5648133,  5.4660729,  7.8242510,  9.2654410,  8.7051126,  9.1776547,  8.7582293, -0.1180924,  3.0975518,  6.2641685,  7.6775858,  9.6217780,  9.9818666,  9.9958351,  9.6264304,  0.1774666,  1.2241896,  4.8679552,  5.2711378,  6.0868630,  9.9227375,  8.1917069,  8.7833352,  8.5495586,  8.7310343), 
                            "score" = c( 0,  3,  5,  8,  9,  9,  9,  9,  0,  3,  6,  8, NA, 10, 10, 10,  0, NA, NA,  4,  6, 10,  8,  9,  9,  9), 
                            "outliers" = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE), 
                            "na.flag" = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE), 
                            "type" = c("signal", "signal", "signal", "signal", "signal", "signal", "signal", "signal", "signal", "signal", "signal", "signal", "signal", "signal", "signal", "signal", "signal", "signal", "signal", "outlier", "signal", "signal", "signal", "signal", "signal", "signal"))

    mod.lmb <- mlmbreak(score ~ 0 + bp(time, "10"), cluster = "ID", data = data.test)
    expect_equal(coef(mod.lmb)$breakpoint, c(69.35013,  67.97003, 101.00000), tol = 1e-3)
})

##----------------------------------------------------------------------
### test-auto-previous-bug.R ends here
