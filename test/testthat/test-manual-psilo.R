### test-manual-psilo.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr  8 2024 (18:40) 
## Version: 
## Last-Updated: apr  8 2024 (20:04) 
##           By: Brice Ozenne
##     Update #: 7
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
    library(ggpubr) 

    library(lmbreak)
}

context("Reproduce previous results")
data("SDIpsilo")


## * graphical display
data(SDIpsilo)
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
if(FALSE){

    SDIpsilo.red <- SDIpsilo[!is.na(SDIpsilo$score) & SDIpsilo$type %in% c("signal","added"),]

    ls.elmbreak <- by(data = SDIpsilo.red, INDICES = SDIpsilo.red$id, FUN = function(iDF){lmbreak(score ~ bp(time, "101"), data = iDF, n.iter = 350)})
    ls.gglmbreak <- lapply(1:length(ls.elmbreak), function(iID){
        autoplot(ls.elmbreak[[iID]], title = paste0("id = ",iID,",cv = ", ls.elmbreak[[iID]]$cv))$plot
    })

    ggall <- do.call(ggarrange,c(ls.gglmbreak, list(common.legend = TRUE, legend = "bottom")))

}

##----------------------------------------------------------------------
### test-manual-psilo.R ends here
