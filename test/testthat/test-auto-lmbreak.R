### test-lmbreak.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Apr  6 2024 (12:19) 
## Version: 
## Last-Updated: apr  8 2024 (17:58) 
##           By: Brice Ozenne
##     Update #: 25
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
    library(segmented)   
    library(lava)   
    library(ggplot2) 

    library(lmbreak)
}

context("Compare to segmented package")

## * Comparison with segmented
## example from the segmented package
set.seed(12)
xx <- 1:100
zz <- runif(100)
yy <- 2+1.5*pmax(xx-35,0)-1.5*pmax(xx-70,0)+15*pmax(zz-.5,0)+rnorm(100,0,2)
dati <- data.frame(x=xx,y=yy,z=zz)
out.lm <- lm(y~x,data=dati)

## ggplot(dati, aes(xx, yy)) + geom_point()

test_that("Pattern 11", {

    ## estimation
    GS.1bp <- segmented(out.lm, npsi = 1)
    test.1bp <- lmbreak(y~bp(x,"11"), data = dati)
    expect_equal(coef(test.1bp, "breakpoint"), GS.1bp$psi[,"Est."], tol = 1e-3)
    expect_equal(test.1bp$breakpoint$Us,"Us1")
    expect_equal(test.1bp$breakpoint$sign,1)

    ## model fit
    expect_equal(as.double(fitted(test.1bp)), as.double(predict(test.1bp, newdata = dati)$estimate), tol = 1e-3)
    df.fit1bp <- predict(test.1bp, newdata = data.frame(x = seq(0,1,length.out =100)), keep.newdata = TRUE, se.fit = FALSE, interval = "confidence")
    plot(test.1bp)
    
})

test_that("Pattern 111", {

    ## estimation
    GS.2bp <- segmented(out.lm, npsi = 2)
    test.2bp <- lmbreak(y~bp(x,"111"), data = dati)                                        
    expect_equal(coef(test.2bp, "breakpoint"), as.double(GS.2bp$psi[,"Est."]), tol = 1e-3)
    expect_equal(test.2bp$breakpoint$Us,c("Us1","Us2"))
    expect_equal(test.2bp$breakpoint$sign,c(1,1))

    ## model fit
    plot(test.2bp)

})

test_that("Pattern 1111", {

    ## estimation
    GS.3bp <- segmented(out.lm, npsi = 3)
    test.3bp <- lmbreak(y~bp(x,"1111"), data = dati)                                        
    ## expect_equal(coef(test.3bp, "breakpoint"), as.double(GS.3bp$psi[,"Est."]), tol = 1e-3)
    expect_equal(test.3bp$breakpoint$Us,c("Us1","Us2","Us3"))
    expect_equal(test.3bp$breakpoint$sign,c(1,1,1))

    ## model fit
    plot(test.3bp)

})

## * with plateau

test_that("Pattern 01", {
    
    set.seed(1)
    df.01 <- simBreak(c(20,50), breakpoint = c(0,1,2), slope = c(0,1), sigma = 0.1)
    test.01 <- lmbreak(Y~bp(X,"01"), data = df.01)
    expect_equal(coef(test.01, type = "breakpoint"), 0.9920744, tol = 1e-3)
    expect_equal(test.01$breakpoint$Us,c("Us1"))
    expect_equal(test.01$breakpoint$sign,1)
    plot(test.01)

    ## convergence issue
    set.seed(1)
    df.01 <- simBreak(c(10,50), breakpoint = c(0,1,2), slope = c(0,1), sigma = 0.1)
    test.01 <- suppressWarnings(lmbreak(Y~bp(X,"01", 1), data = df.01))
    plot(test.01)

    GS.1bp <- segmented(lm(Y~X, data = df.01), npsi = 1)
    test.11 <- lmbreak(Y~bp(X,"11"), data = df.01)

    expect_equal(coef(test.01, type = "breakpoint"), 1.001113, tol = 1e-3)
    expect_equal(coef(test.11, type = "breakpoint"), GS.1bp$psi[,"Est."], tol = 1e-3)

})

test_that("Pattern 10", {
    
    set.seed(1)
    df.10 <- simBreak(c(20,50), breakpoint = c(0,1,2), slope = c(1,0), sigma = 0.1)
    test.10 <- lmbreak(Y~bp(X,"10"), data = df.10)
    expect_equal(test.01$breakpoint$Us,c("Us1"))
    expect_equal(test.01$breakpoint$sign,1)

    plot(test.10)
    expect_equal(coef(test.10, type = "breakpoint"), 0.9881715, tol = 1e-3)

})

test_that("Pattern 010", {

    ## estimation
    test.010 <- lmbreak(y~bp(x,"010"), data = dati)
    expect_equal(coef(test.010, "breakpoint"), c(33.78736,70.10277), tol = 1e-3)
    expect_equal(test.010$breakpoint$Us,c("I(Us1 - Us2)","I(Us1 - Us2)"))
    expect_equal(test.010$breakpoint$sign,c(1,-1))

    ## model fit
    plot(test.010)
})

test_that("Pattern 011", {

    ## estimation
    test.011 <- lmbreak(y~bp(x,"011"), data = dati)
    expect_equal(coef(test.011, "breakpoint"), c(33.67924, 71.77749), tol = 1e-3)
    expect_equal(test.011$breakpoint$Us,c("Us1","Us2"))
    expect_equal(test.011$breakpoint$sign,c(1,1))

    ## model fit
    plot(test.011)

})

test_that("Pattern 011", {

    ## estimation
    test.110 <- suppressWarnings(lmbreak(y~bp(x,"110"), data = dati))
    expect_equal(coef(test.110, "breakpoint"), c(32.68963358, 70.27439266), tol = 1e-3)
    expect_equal(test.110$breakpoint$Us,c("I(Us1 - Us2)","I(Us1 - Us2)"))
    expect_equal(test.110$breakpoint$sign,c(1,-1))

    ## model fit
    plot(test.110)

})

test_that("Pattern 101", {
    
    set.seed(1)
    df.101 <- simBreak(c(20,50), breakpoint = c(0,1,2,3), slope = c(1,0,-1), sigma = 0.1)
    test.101 <- lmbreak(Y~bp(X,"101"), data = df.101)
    expect_equal(coef(test.101, type = "breakpoint"), c(0.9898938, 1.9811977), tol = 1e-3)

    expect_equal(test.101$breakpoint$Us,c("I(Us0 - Us1)","Us2"))
    expect_equal(test.101$breakpoint$sign,c(-1,1))

    plot(test.101)

    set.seed(1)
    df.101 <- simBreak(c(10,50), breakpoint = c(0,1,2,4), slope = c(2,0,-1), sigma = 0.1)
    test.101 <- suppressWarnings(lmbreak(Y~bp(X,"101"), data = df.101))
    expect_equal(coef(test.101, type = "breakpoint"), c(1.004624, 1.991102), tol = 1e-3)
    plot(test.101)
})


test_that("Pattern 1010", {
    
    set.seed(10)
    df.1010 <- simBreak(c(20,50), breakpoint = c(0,1,2,3,4), slope = c(1,0,-0.9,0), sigma = 0.05)
    test.1010 <- lmbreak(Y~bp(X,"1010"), data = df.1010)
    expect_equal(coef(test.1010, type = "breakpoint"), c(1.009588, 1.991627, 3.011256), tol = 1e-3)

    plot(test.1010)

    set.seed(10)
    df.1010 <- simBreak(c(20,50), breakpoint = c(0,1,2,3,4), slope = c(1,0,-0.9,0), sigma = 0.1)
    test.1010 <- lmbreak(Y~bp(X,"1010", init = 1:3), data = df.1010)
    expect_equal(coef(test.1010, type = "breakpoint"), c(1.019542, 1.982325, 3.021334), tol = 1e-3)

})
##----------------------------------------------------------------------
### test-lmbreak.R ends here
