# Load required packages

library(readxl)
source("include.R")


# Read in the data and give handy variable names

cases <- as.data.frame(read_excel("input/Line listing_CY.xlsx"))
vacc <- as.data.frame(read_excel("input/VACCINATION DATA UNTIL 17 03 2022.xlsx"))
pop <- as.data.frame(read_excel("input/Population 2019.xlsx", skip=1))
pop <- pop[-c(1,nrow(pop)), ] 
names(cases) <- c("oid", "age", "dres", "dsamp", "hosp", "ddeath", "death", "cdeath", "vaccdate1", "vaccdate2", "vaccdate3", "vaccType.raw")
names(vacc) <- c("date", "vaccType", "dose", "n", "agegrp")
names(pop) <- c("agegrp", "n")


# Define correspondence with vaccine codes & other categorical variables
vaccMat <- c("PFIZER, BIONTECH" = "Pfizer", "COMIRNATY" = "Pfizer", 
  "ASTRAZENECA, UNIV. OF OXFORD" = "AZ", "MODERNA" = "Moderna", "JOHNSON & JOHNSON" = "J&J")
vaccMat2 <- c("ASTRAZEN01" = "AZ", "COMIRNATY01" = "Pfizer", "JOHNSON01" = "J&J", "MODERNA01" = "Moderna", "PFIZER01" = "Pfizer")
vaccTypes <- c("Pfizer", "Moderna", "AZ", "J&J", "other")
vacclbl <- function(i=NULL) {
  a <- c("Unvaccinated", "1-dose vaccinated", "2-dose vaccinated", "3-dose vaccinated")
  if (is.null(i)) return(a)
  a[i+1]
}


# Housekeeping and recoding variables

# cases
cases$hosp <- cases$hosp == sort(unique(cases$hosp))[1]
cases$death <- cases$death == sort(unique(cases$death))[1]
cases$dres <- as.Date(cases$dres)
cases$dsamp <- as.Date(cases$dsamp)
cases$ddeath <- as.Date(cases$ddeath)
cases$vaccdate1 <- as.Date(cases$vaccdate1)
cases$vaccdate2 <- as.Date(cases$vaccdate2)
cases$vaccdate3 <- as.Date(cases$vaccdate3)
cases$vaccType <- unname(vaccMat[cases$vaccType.raw])
cases$vaccType[is.na(cases$vaccType) & !is.na(cases$vaccType.raw)] <- "other"
cases$vaccType <- factor(cases$vaccType, levels=vaccTypes)
# vacc
vacc$date <- as.Date(vacc$date)
vacc$vaccType <- unname(vaccMat2[vacc$vaccType])
vacc$vaccType[is.na(vacc$vaccType)] <- "other"
vacc$vaccType <- factor(vacc$vaccType, levels=vaccTypes)
vacc$agegrp <- factor(vacc$agegrp, levels=sort(pop$agegrp))
# pop
pop$agegrp <- factor(pop$agegrp, levels=sort(pop$agegrp))
rownames(pop) <- NULL


# Handle likely errors and non-sensical values

cases$dsamp[with(cases, as.integer(dres-dsamp)==28)] <- cases$dsamp[with(cases, as.integer(dres-dsamp)==28)] + 28
cases$dsamp[with(cases, as.integer(dres-dsamp)==30)] <- cases$dsamp[with(cases, as.integer(dres-dsamp)==30)] + 30
cases$dsamp[with(cases, as.integer(dres-dsamp)==31)] <- cases$dsamp[with(cases, as.integer(dres-dsamp)==31)] + 31
cases$dsamp[with(cases, as.integer(dres-dsamp)==365)] <- cases$dsamp[with(cases, as.integer(dres-dsamp)==365)] + 365
cases$date <- cases$dsamp
cases$date[with(cases, as.integer(dres-dsamp)>31)] <- cases$dres[with(cases, as.integer(dres-dsamp)>31)]
vacc <- subset(vacc, date>="2020-12-27")
vacc[with(vacc, which(dose==3 & date==min(date[dose==3]))),"dose"] <- 1   # Erroneous 3rd dose at 20/1/2020, probably is 1st...

# Ugly TEMPORARY population hack to make the numbers have sense. TO BE REMOVED ONCE CORRECT VACC DATA RECEIVED
pop$n[13:17] <- round(pop$n[13:17] * c(1.05, 1.12, 1.12, 1.2, 1.3))


# Exclude cases vaccinated with other vaccines than the licensed ones.
out <- list()  # Placeholder for various output used in the report
out$exclCasesUnkVacc <- subset(cases, vaccType=="other")  # Keep the number excluded for reference
cases <- subset(cases, is.na(vaccType) | vaccType!="other")


# Calculate derivative variables

cases$agegrp <- cut(cases$age, breaks=c(5,12,18,seq(25,90,5),200), right=FALSE, labels=levels(pop$age))
cases$vaccStatus <- factor(vacclbl(0), levels=vacclbl())
# Vaccinated status is given 7 days after administration of the respective dose
cases$vaccStatus[!is.na(cases$vaccdate1) & cases$vaccdate1<=(cases$date-7)] <- vacclbl(1)
cases$vaccStatus[!is.na(cases$vaccdate2) & cases$vaccdate2<=(cases$date-7)] <- vacclbl(2)
cases$vaccStatus[!is.na(cases$vaccdate3) & cases$vaccdate3<=(cases$date-7)] <- vacclbl(3)
cases$tdiff <- 0   # Delay between vaccination and becoming a case
cases$tdiff[cases$vaccStatus==vacclbl(1)] <- with(cases, as.integer(date[vaccStatus==vacclbl(1)] - vaccdate1[vaccStatus==vacclbl(1)] - 7))
cases$tdiff[cases$vaccStatus==vacclbl(2)] <- with(cases, as.integer(date[vaccStatus==vacclbl(2)] - vaccdate2[vaccStatus==vacclbl(2)] - 7))
cases$tdiff[cases$vaccStatus==vacclbl(3)] <- with(cases, as.integer(date[vaccStatus==vacclbl(3)] - vaccdate3[vaccStatus==vacclbl(3)] - 7))


# Aggregate vaccination data and case data

# First aggregate doses administered in each date and group
vaccAggr1 <- aggregate(vacc[,"n",drop=FALSE], vacc[,c("dose","date","agegrp")], sum)
vaccAggr1 <- reshape(vaccAggr1, idvar=c("date","agegrp"), timevar="dose", direction="wide")
vaccAggr1 <- merge(expand.grid(date=seq.Date(min(vacc$date), max(vacc$date), by="day"), agegrp=levels(vacc$agegrp)), vaccAggr1, all=TRUE)  # Ensure all dates are included
vaccAggr1$n.1[is.na(vaccAggr1$n.1)] <- 0
vaccAggr1$n.2[is.na(vaccAggr1$n.2)] <- 0
vaccAggr1$n.3[is.na(vaccAggr1$n.3)] <- 0
vaccAggr1$unvacc <- pop$n[match(vaccAggr1$agegrp, pop$agegrp)]

# Then calculate the "pool" of people in each vaccination status, per age group
vaccAggr2 <- do.call(rbind, by(vaccAggr1, vaccAggr1[,"agegrp"], function(x) {
  x$n.1 <- vecshift(cumsum(x$n.1), -7)
  x$n.2 <- vecshift(cumsum(x$n.2), -7)
  x$n.3 <- vecshift(cumsum(x$n.3), -7)
  x$vacc3 <- x$n.3
  x$vacc2 <- x$n.2 - x$vacc3
  x$vacc1 <- x$n.1 - x$vacc2 - x$vacc3
  x$unvacc <- x$unvacc - x$vacc1 - x$vacc2 - x$vacc3
  x <- x[,c("date","agegrp","unvacc","vacc1","vacc2","vacc3")]
  x <- reshape(x, idvar=c("date","agegrp"), varying=c("unvacc","vacc1","vacc2","vacc3"), v.names="n", times=vacclbl(), timevar="vaccStatus", direction="long")
  x$vaccStatus <- factor(x$vaccStatus, levels=levels(cases$vaccStatus))
  x
}))
rownames(vaccAggr2) <- NULL

# Aggregate case data, and merge with vaccination data
cases$case <- 1
casesAggr <- aggregate(cases[,c("case","hosp","death")], cases[,c("date","agegrp","vaccStatus")], sum)
dat <- merge(vaccAggr2, casesAggr, all.x=TRUE)
dat$case[is.na(dat$case)] <- 0
dat$hosp[is.na(dat$hosp)] <- 0
dat$death[is.na(dat$death)] <- 0
dat$agegrp <- relevel(dat$agegrp, "50-54")
dat$isoweek <- isoweek(dat$date)
dat <- subset(dat, n>0)  # Exclude dates where there are no people in a vaccinated "pool"

# Aggregate also across all age groups (to be able to obtain non-age-adjusted estimates)
datU <- aggregate(dat[,c("n","case","hosp","death")], dat[,c("date","vaccStatus","isoweek")], sum)
datU <- datU[order(datU$date),]

