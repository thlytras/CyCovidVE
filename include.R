# Make nice color palettes available
library(RColorBrewer)


# Functions to make different functionalities available


# Function to shift a vector elements d places
vecshift <- function(x, d=0) {
  if (d==0) return(x)
  if (d>0) return(c(x[-(1:d)], rep(0,d)))
  if (d<0) return(c(rep(0, -d), x[1:(length(x)+d)]))
}


# Function to calculate ISO weeks
isoweek <- function(x, type="both_num", sep="-", inv=FALSE, colnames=c("isoyear","isoweek")) {
  alts=c("week","year","both_text","both_num","matrix")
  if(!(type %in% alts)) stop("Unknown isoweek type requested!")
  x.date<-as.Date(x)
  x.weekday<-as.integer(format(x.date,"%w"))
  x.weekday[x.weekday==0]=7
  x.nearest.thu<-x.date-x.weekday+4
  x.isoyear<-as.integer(substring(x.nearest.thu,1,4)) # Μπορεί οι πρώτες μέρες του χρόνου να ανήκουν (κατά ISO) στην προηγούμενη χρονιά!
  x.isoweek<-(as.integer(x.nearest.thu-as.Date(paste(x.isoyear,"-1-1",sep="")))%/%7)+1
  switch(type,
    week = x.isoweek,
    year = x.isoyear,
    both_text = if (inv) {
      ifelse((is.na(x.isoyear) | is.na(x.isoweek)),NA,paste(x.isoweek,x.isoyear,sep=sep))
    } else {
      ifelse((is.na(x.isoyear) | is.na(x.isoweek)),NA,paste(x.isoyear,x.isoweek,sep=sep))
    },
    both_num = ifelse((is.na(x.isoyear) | is.na(x.isoweek)),NA,x.isoyear*100+x.isoweek),
    matrix = if (inv) {
      `colnames<-`(cbind(x.isoweek, x.isoyear), rev(colnames))
    } else {
      `colnames<-`(cbind(x.isoyear, x.isoweek), colnames)
    }
  )
}


# Function to calculate age-adjusted Kaplan Meier estimates
ageAdjKM <- function(v, outc, from=min(dat$date), to=max(dat$date)) {
  n <- lapply(levels(dat$agegrp), function(a) subset(dat, vaccStatus==v & agegrp==a & date>=from & date<=to)[,c("date", "n")])
  n <- suppressWarnings(Reduce(function(x,y) merge(x,y,all=TRUE, by="date"), n))
  colnames(n)[-1] <- paste0("n", levels(dat$agegrp))
  x <- lapply(levels(dat$agegrp), function(a) subset(dat, vaccStatus==v & agegrp==a & date>=from & date<=to)[,c("date", outc)])
  x <- suppressWarnings(Reduce(function(x,y) merge(x,y,all=TRUE, by="date"), x))
  colnames(x)[-1] <- paste0("x", levels(dat$agegrp))
  dates <- n$date
  n <- n[,-1]; x <- x[,-1]
  pp <- matrix(pop$n[match(levels(dat$agegrp), pop$agegrp)], ncol=length(levels(dat$agegrp)), nrow=nrow(n), byrow=TRUE)
  pp[is.na(x)] <- NA
  w <- pp/rowSums(pp, na.rm=TRUE)
  qstar <- rowSums((x/n)*w, na.rm=TRUE)
  varlogS <- cumsum(rowSums(w^2 * (1-x/n)*(x/n)/n, na.rm=TRUE) * (1/(1-qstar)^2))
  res <- data.frame(date = dates, S = cumprod(1-qstar), se.logS = sqrt(varlogS))
  res$S.lo <- with(res, exp(log(S)-1.96*se.logS))
  res$S.hi <- with(res, exp(log(S)+1.96*se.logS))
  res$S.lo[res$S.lo>1] <- 1
  res
}


# Function to calculate simple Kaplan-Meier estimates (non-age-adjusted)
simpleKM <- function(v, outc, from=min(datU$date), to=max(dat$date)) {
  x <- subset(datU, vaccStatus==v & date>=from & date<=to)[,c("date","n",outc)]
  x$S <- cumprod(1-x[,outc]/x$n)
  x$varS <- with(x, S^2*cumsum(get(outc)/(n*(n-get(outc)))))
  x$varlogS <- x$varS / x$S^2
  x$S.lo <- with(x, S*exp(1.96*sqrt(varlogS)))
  x$S.hi <- with(x, S*exp(-1.96*sqrt(varlogS)))
  x$S.lo[x$S.lo>1] <- 1
  return(x)
}


# Function to add alpha transparency to a color (I like semitransparent bands in figures!)
addalpha <- function(colors, alpha=1.0) {
  r <- col2rgb(colors, alpha=T)
  # Apply alpha
  r[4,] <- alpha*255
  r <- r/255.0
  return(rgb(r[1,], r[2,], r[3,], r[4,]))
}


# Function to draw a polygon (band) with Kaplan-Meier-like "steps"
# (e.g. to indicate the confidence bands of the Kaplan-Meier estimator)
polygonS <- function(x, lo, hi, col) {
  polygon(
    c(rep(x, each=2)[-c(1)], rev(rep(x, each=2)[-c(1)])),
    c(rep(lo, each=2)[-c(length(lo)*2)], rev(rep(hi, each=2)[-c(length(hi)*2)])),
    border=NA, col=col
  )
}


# Make a descriptive table of rates
makeRatesTable <- function(from, to) {
  # Function to estimate the median age in a distribution of agegroups, using a loess fit
  findMedianAge <- function(x) {
    x[is.na(x)] <- 0
    midp <- c(8, 14.5, 21, seq(27, 87, 5), 92)
    p <- predict(loess(x~midp), 1:100)
    p[is.na(p) | p<0] <- 0
    which(cumsum(p)/sum(p,na.rm=TRUE)>0.5)[1]
  }
  pr <- function(x, f="%.1f") { # Simple pretty-printer
    if (is.na(x) | x==0) return("-")
    sprintf(f, x)
  }
  sapply(levels(cases$vaccStatus), function(v) {
    s <- subset(cases, date>=from & date<=to & vaccStatus==v)
    c.a <- tapply(s$case, s$agegrp, sum)
    c.a[is.na(c.a)] <- 0
    h.a <- tapply(s$hosp, s$agegrp, sum)
    h.a[is.na(h.a)] <- 0
    d.a <- tapply(s$death, s$agegrp, sum)
    d.a[is.na(d.a)] <- 0
    n.a <- with(subset(vaccAggr2, date>=from & date<=to & vaccStatus==v), tapply(n, agegrp, sum))
    n.a[n.a==0] <- NA
    res <- c(
      "pop" = round(sum(n.a, na.rm=TRUE)/as.integer(as.Date(to)-as.Date(from))),
      "mdAge" = if (round(sum(n.a, na.rm=TRUE)/as.integer(as.Date(to)-as.Date(from)))>0) findMedianAge(n.a) else "-",
      "cases" = sum(c.a),
      "caseR" = pr(sum(c.a)/(sum(n.a, na.rm=TRUE)/365.25*12)*100000),
      "caseAR" = pr(sum((c.a/(n.a/365.25*12)*100000)*pop$n/sum(pop$n), na.rm=TRUE)),
      "hosp" = sum(h.a),
      "hospR" = pr(sum(h.a)/(sum(n.a, na.rm=TRUE)/365.25*12)*100000),
      "hospAR" = pr(sum((h.a/(n.a/365.25*12)*100000)*pop$n/sum(pop$n), na.rm=TRUE)),
      "death" = sum(d.a),
      "deathR" = pr(sum(d.a)/(sum(n.a, na.rm=TRUE)/365.25*12)*100000),
      "deathAR" = pr(sum((d.a/(n.a/365.25*12)*100000)*pop$n/sum(pop$n), na.rm=TRUE))
    )
  })
}


# Function to plot a set of Kaplan-Meier curves (for any outcome, any period, etc)
plotAllKM <- function(outc="case", from="2021-12-20", to=max(dat$date), show=1:3, adj=TRUE, inc=TRUE, alpha=0.1, title="") {
  cols <- brewer.pal(4, "Spectral")    # cols <- c("red", "orange", "green", "blue")
  vcs <- c("Unvaccinated", paste0(1:3, "-dose vaccinated"))
  show <- c(0, show[show %in% c(1:3)])
  if (adj) {
    with(ageAdjKM("Unvaccinated", outc, from=from, to=to), 
      plot(date, if (inc) (1-S) else S, type="n", col="red", xlab="Ημερομηνία", yaxt="n", bty="l", main=title,
        ylab=if(inc) "Αθροιστική επίπτωση (%)" else "Ποσοστό επιβίωσης (%)"))
    axis(2, at=axTicks(2), labels=paste0(axTicks(2)*100, ""), las=2)
    sapply(show, function(i) with(ageAdjKM(vcs[i+1], outc, from=from, to=to), {
      polygonS(date, if (inc) (1-S.hi) else S.hi, if (inc) (1-S.lo) else S.lo, col=addalpha(cols[i+1], alpha))
      points(date, if (inc) (1-S) else S, type="s", col=cols[i+1])
    }))
  } else {    
    with(simpleKM("Unvaccinated", outc, from=from, to=to), 
      plot(date, if (inc) (1-S) else S, type="n", col="red", xlab="Ημερομηνία", yaxt="n", bty="l", main=title,
        ylab=if(inc) "Αθροιστική επίπτωση (%)" else "Ποσοστό επιβίωσης (%)"))
    axis(2, at=axTicks(2), labels=paste0(axTicks(2)*100, ""), las=2)
    sapply(show, function(i) with(simpleKM(vcs[i+1], outc, from=from, to=to), {
      polygonS(date, if (inc) (1-S.hi) else S.hi, if (inc) (1-S.lo) else S.lo, col=addalpha(cols[i+1], alpha))
      points(date, if (inc) (1-S) else S, type="s", col=cols[i+1])
    }))
  }
  legend(if (inc) "topleft" else "bottomleft",
    c("Ανεμβολίαστοι", paste("Εμβολιασμένοι με", c("1 δόση", "2 δόσεις", "3 δόσεις")))[show+1], 
    col=cols[show+1], bty="n", lwd=1, seg.len=4, cex=0.8)
}


# Extract age group RRs from a multivariate Poisson model
getRRage <- function(m) {
  res <- exp(coef(summary(m))[grep("^agegrp", names(coef(m))), 1:2] %*% rbind(1, c(0, -1.96, 1.96)))
  colnames(res) <- c("RR", "lo", "hi")
  rownames(res) <- gsub("^agegrp", "", rownames(res))
  res
}


# Draw nice 95% CI bars in a plot (to illustrate VE)
plotBars <- function(x, lo, hi, y=NULL, lwd=1, cex=1, log=TRUE) {
  if (log) log <- "x" else log <- ""
  if (is.null(y)) {
    y <- rev(1:length(x))
  }
  for (i in 1:length(y)) lines(c(lo[i], hi[i]), rep(y[i],2), lwd=lwd)
  points(x, y, pch=15, cex=cex)
}


# Plot RRs per age group
plotAgeAdjRisk <- function(m, log=TRUE, mar=c(5,5,3,9), reflbl="Κατηγορία αναφοράς") {
  x <- getRRage(m)
  i <- which(rownames(x)=="45-49")
  x <- as.data.frame.matrix(rbind(x[1:i,], 1, x[(i+1):nrow(x),]))
  Min <- min(x$lo)
  if (Min>0.01) Min <- floor(Min*100)/100
  if (Min>0.1) Min <- floor(Min*10)/10
  rownames(x)[i+1] <- "50-54"
  par(family="Fira Sans", mar=mar)
  plot(0, type="n", bty="n", axes=FALSE, xlab="Σχετικός κίνδυνος (Rate Ratio)", ylim=c(0.5,nrow(x)), xlim=c(Min, max(x$hi)+0.2+(max(x$hi)>5)*2), log=c("","x")[log+1], ylab=NA)
  abline(v=1, lty="dotted")
  plotBars(x$RR, x$lo, x$hi)
  axis(1, at=axTicks(1)[seq(1,length(axTicks(1)),by=2)], labels=sprintf("%s",axTicks(1))[seq(1,length(axTicks(1)),by=2)])
  if ((length(axTicks(1)))>=4) axis(1, at=axTicks(1)[seq(4,length(axTicks(1)),by=2)], labels=sprintf("%s",axTicks(1))[seq(4,length(axTicks(1)),by=2)])
#  axis(1, at=c(0.01,0.02,0.05,0.1,0.2,0.5,1,2,5,10,20,50), labels=sprintf("%s",c(0.01,0.02,0.05,0.1,0.2,0.5,1,2,5,10,20,50)))
  mtext(c("Ηλικία", rownames(x)), at=(nrow(x)+1):1, side=2, las=1, adj=0.5, line=2, xpd=TRUE)
  rrlab <- c(with(x[1:i,], sprintf(ifelse(hi<10, "%.2f (%.2f–%.2f)", "%.1f (%.1f–%.1f)"), RR, lo, hi)), 
    reflbl, with(x[(i+2):nrow(x),], sprintf(ifelse(hi<10, "%.2f (%.2f–%.2f)", "%.1f (%.1f–%.1f)"), RR, lo, hi)))
  mtext(c("RR", rrlab), at=(nrow(x)+1):1, side=4, adj=0, las=1)
}


# Get predicted number of events from quasi-Poisson model assuming no vaccination
getEventsU <- function(m) {
  NA0 <- function(x) { x[is.na(x)] <- 0; x }
  mfU <- model.matrix(m)
  mfU[,grep("vaccStatus",colnames(mfU))] <- 0
  data.frame(
    x = c(exp((mfU %*% NA0(coef(m)[colnames(mfU)])) + m$offset)),
    n = exp(m$offset),
    var = diag(mfU %*% tcrossprod(NA0(vcov(m)[match(colnames(mfU), rownames(vcov(m))), match(colnames(mfU), colnames(vcov(m)))]), mfU)),
    o = fitted(m)
  )
  # For the variance of the sum: sum(var*x^2)
  # Overdispersion: max(1,sum(m$weights * m$residuals^2)/m$df.r)
}

