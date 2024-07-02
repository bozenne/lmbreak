### test-auto-mlmbreak.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr  8 2024 (18:40) 
## Version: 
## Last-Updated: jul  2 2024 (14:30) 
##           By: Brice Ozenne
##     Update #: 29
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
    library(ggplot2) 

    library(lmbreak)
}

context("mlmbreak on real data in presence of NA")
data("SDIpsilo", package = "lmbreak")


## * graphical display
if(FALSE){

    ggTraj <- ggplot(mapping = aes(x = time,y = score, group = id))
    ggTraj <- ggTraj + geom_line(data = SDIpsilo[!is.na(SDIpsilo$score) & SDIpsilo$type %in% c("signal","added"),], color = palette.colors()["skyblue"])
    ggTraj <- ggTraj + geom_point(data = SDIpsilo[SDIpsilo$type == "signal",], color = palette.colors()["skyblue"], shape = 19)
    ggTraj <- ggTraj + geom_point(data = SDIpsilo[SDIpsilo$type == "noise",], color = palette.colors()["orange"], shape = 8)
    ggTraj <- ggTraj + geom_point(data = SDIpsilo[SDIpsilo$type == "added",], color = palette.colors()["vermillion"], shape = 17)
    ggTraj <- ggTraj + geom_point(data = SDIpsilo[SDIpsilo$type == "trailing",], color = palette.colors()["gray"], shape = 15)
    ggTraj <- ggTraj + facet_wrap(~id, labeller = label_both)
    ggTraj <- ggTraj + coord_cartesian(ylim = c(-1.25,10))
    ggTraj <- ggTraj + xlab("time") + guides(color = "none")
    ggTraj

}

## * lmbreak
test_that("mlmbreak with NA", {

    SDIpsilo.red <- SDIpsilo[SDIpsilo$type %in% c("signal","added"),]
    e.mlmbreak <- mlmbreak(score ~ 0 + bp(time, "101"), cluster = "id", data = SDIpsilo.red, trace = FALSE)
    expect_equal(logLik(e.mlmbreak), -266.5181, tol = 1e-3)
    suppressWarnings(plot(e.mlmbreak, ylim = c(0,10)))

    SDIpsilo.redNNA <- SDIpsilo.red[!is.na(SDIpsilo.red$score),]
    eNNA.mlmbreak <- mlmbreak(score ~ 0 + bp(time, "101"), cluster = "id", data = SDIpsilo.redNNA, trace = FALSE)

    test <- model.tables(e.mlmbreak, format = "array")
    GS <- model.tables(eNNA.mlmbreak, format = "array")
    expect_equal(test,GS, tol = 1e-3)
    expect_equal(logLik(e.mlmbreak), logLik(eNNA.mlmbreak), tol = 1e-3)
    
})

##----------------------------------------------------------------------
### test-auto-mlmbreak.R ends here
