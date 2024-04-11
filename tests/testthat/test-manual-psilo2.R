### test-manual-psilo2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr  8 2024 (20:05) 
## Version: 
## Last-Updated: apr 11 2024 (10:36) 
##           By: Brice Ozenne
##     Update #: 12
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:


if(FALSE){
    library(readxl)
    library(data.table)
    library(ggplot2)

    ## * load data
    dtW <- as.data.table(read_excel("SDI_Psilo_22.11.2023.xlsx"))
    dtL <- melt(dtW, id.vars = c("Study","PatientID"), variable.name = "time", value.name = "score")
    dtL$time.num <- as.numeric(as.character(dtL$time))

    ## * display
    ggTraj <- ggplot(dtL, aes(x = time.num, y = score, group = PatientID, color = Study))
    ggTraj <- ggTraj + geom_line() + geom_point()
    ggTraj <- ggTraj + facet_wrap(~PatientID)
    ggTraj

    ## * Analyse only SERT individuals
    dtL.SERT <- dtL[grep("SERT",dtL$PatientID),,drop=FALSE]
    ggTraj %+% dtL.SERT

    e.mlmbreak1010 <- mlmbreak(score ~ 0 + bp(time.num, c("1010","101","11")), cluster = "PatientID", data = dtL.SERT)
    e.mlmbreak1010
    plot(e.mlmbreak1010, labeller = label_value, scales = "free", ylim = c(0,100))
    model.tables(e.mlmbreak1010, format = "list")

    library(segmented)
    ls.lm <- dtL.SERT[,.(.(lm(score ~ time.num, data = .SD))), by = "PatientID"]$V1

    segmented(ls.lm[[1]], seg.Z=~time.num, npsi = 2)
}


##----------------------------------------------------------------------
### test-manual-psilo2.R ends here
