### doc-data.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr  8 2024 (18:50) 
## Version: 
## Last-Updated: apr  8 2024 (19:50) 
##           By: Brice Ozenne
##     Update #: 14
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * SDIpsilo
## SDIpsilo <- as.data.frame(dtL.data[,.(id = as.numeric(as.factor(CIMBI.ID)),time=time.num,type=as.character(datapoint),score)][id %in% c(5,7,17)==FALSE])
## SDIpsilo$id <- as.factor(as.numeric(as.factor(SDIpsilo$id)))
## SDIpsilo <- rbind(SDIpsilo, data.frame(id = 3, time = 0, type = "added", score = 0))
## SDIpsilo <- SDIpsilo[order(SDIpsilo$id,SDIpsilo$time),]
## SDIpsilo[SDIpsilo$id == 1 & SDIpsilo$time == 140,"type"] <- "noise"
## SDIpsilo[SDIpsilo$id == 2 & SDIpsilo$time %in% c(280,300),"type"] <- "trailing"
## SDIpsilo[SDIpsilo$id == 3 & SDIpsilo$time %in% c(300),"type"] <- "trailing"
## SDIpsilo[SDIpsilo$id == 4 & SDIpsilo$time %in% c(440),"type"] <- "trailing"
## SDIpsilo[SDIpsilo$id == 6 & SDIpsilo$time %in% c(260),"type"] <- "trailing"
## SDIpsilo[SDIpsilo$id == 7 & SDIpsilo$time %in% c(320,340),"type"] <- "trailing"
## SDIpsilo[SDIpsilo$id == 8 & SDIpsilo$time %in% c(300,320),"type"] <- "trailing"
## SDIpsilo[SDIpsilo$id == 9 & SDIpsilo$time %in% c(140),"type"] <- "noise"
## SDIpsilo[SDIpsilo$id == 10 & SDIpsilo$time %in% c(340),"type"] <- "trailing"
## SDIpsilo[SDIpsilo$id == 11 & SDIpsilo$time %in% c(320),"type"] <- "trailing"
## SDIpsilo[SDIpsilo$id == 13 & SDIpsilo$time %in% c(320,340,360),"type"] <- "trailing"
## SDIpsilo[SDIpsilo$id == 13 & SDIpsilo$time %in% c(180),"type"] <- "noise"
## SDIpsilo[SDIpsilo$id == 14 & SDIpsilo$time %in% c(20,40),"type"] <- "trailing"
## SDIpsilo[SDIpsilo$id == 15 & SDIpsilo$time %in% c(220),"type"] <- "noise"
## 
## save(SDIpsilo, file = "data/SDIpsilo.rda")
#' @title SDI data after Psilocibin
#' @name SDIpsilo
#' @rdname data-SDIpsilo
#'
#' @description Data from an experimental study where 16 healthy subjects receve an oral dose of psilocibin and their SDI (subjective drug intensity) score is followed over time.
#'
#' \itemize{
#' \item \code{id}: study participant.
#' \item \code{time}: time (in minutes) since psilocibin intake.
#' \item \code{datapoint}: \code{"signal"} or \code{"noise"}. Indicates when the SDI measurement was done just after the participant visited the toilet which substantively, but temporarily, disrupted the experience (noise)
#' \item \code{score}: SDI score between 0 and 10
#' }
#' 
#' @docType data
#' @usage data(SDIpsilo)
#' @references Stenbæk DS, Madsen MK, Ozenne B, et al. Brain serotonin 2A receptor binding predicts subjective temporal and mystical effects of psilocybin in healthy humans. Journal of Psychopharmacology. 2021;35(4):459-468. doi:10.1177/0269881120959609.
#' @keywords datasets
NULL



##----------------------------------------------------------------------
### doc-data.R ends here
