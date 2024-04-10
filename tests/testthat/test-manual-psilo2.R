### test-manual-psilo2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr  8 2024 (20:05) 
## Version: 
## Last-Updated: apr 10 2024 (15:50) 
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


library(readxl)
library(data.table)
library(ggplot2)

## * load data
if(FALSE){
    dtW <- as.data.table(read_excel("SDI_Psilo_22.11.2023.xlsx"))
    dtL <- melt(dtW, id.vars = c("Study","PatientID"), variable.name = "time", value.name = "score")
    dtL$time.num <- as.numeric(as.character(dtL$time))
}

## * display
if(FALSE){
    ggTraj <- ggplot(dtL, aes(x = time.num, y = score, group = PatientID, color = Study))
    ggTraj <- ggTraj + geom_line() + geom_point()
    ggTraj <- ggTraj + facet_wrap(~PatientID)
    ggTraj
}

## * Analyse only SERT individuals
if(FALSE){
    dtL.SERT <- dtL[grep("SERT",dtL$PatientID),,drop=FALSE]
    ggTraj %+% dtL.SERT

    e.mlmbreak1010 <- mlmbreak(score ~ 0 + bp(time.num, c("1010","101","11")), cluster = "PatientID", data = dtL.SERT)
    e.mlmbreak1010
    plot(e.mlmbreak1010, labeller = label_value, scales = "free", ylim = c(0,100))
    model.tables(e.mlmbreak1010, format = "array")
    
    ## e21.lmbreak1010 <- lmbreak(score ~ 0 + bp(time.num, c("1010"), c(150,175,210)), data = dtL.SERT[!is.na(dtL.SERT$score) & dtL.SERT$PatientID=="SERT-Psilo-021",])
    ## plot(e21.lmbreak1010)
}


##----------------------------------------------------------------------
### test-manual-psilo2.R ends here
