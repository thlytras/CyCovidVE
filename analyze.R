# Load required packages
library(plotrix)

# If calling without having called prep.R first, load the saved prepared data
if (!exists("dat") && !exists("out")) load("intermediate.RData")



# Fit multivariate models
cat("Fitting models...\n")
m <- list(
  death = glm(death ~ vaccStatus + agegrp + factor(isoweek), offset=log(n), data=subset(dat, agegrp!="05-11" & agegrp!="18-24"), family="quasipoisson"),
  hosp = glm(hosp ~ vaccStatus + agegrp + factor(isoweek), offset=log(n), data=dat, family="quasipoisson"),
  case = glm(case ~ vaccStatus + agegrp + factor(isoweek), offset=log(n), data=dat, family="quasipoisson")
)
# Obtaining predicted values for "fully unvaccinated"
cat("Obtaining predicted values for \"fully unvaccinated\", please wait...\n")
predm <- lapply(m, getEventsU)


cat("Calculating all the rest...\n")

# Make descriptive rate tables
out$rtab1 <- makeRatesTable("2020-12-27", "2021-6-20")
out$rtab2 <- makeRatesTable("2021-6-21", "2021-12-19")
out$rtab3 <- makeRatesTable("2021-12-20", max(dat$date))
colnames(out$rtab1) <- colnames(out$rtab2) <- colnames(out$rtab3) <- 
  c("Ανεμβολίαστοι", "Εμβολιασμένοι με 1 δόση", "Εμβολιασμένοι με 2 δόσεις", "Εμβολιασμένοι με 3 δόσεις")
rownames(out$rtab1) <- rownames(out$rtab2) <- rownames(out$rtab3) <- 
  c("Πληθυσμός*", "Διάμεση ηλικία", 
  "Κρούσματα", "Κρούσματα / 1000 πληθυσμού", "Κρούσματα / 1000 πληθυσμού, με ηλικιακή στάθμιση",
  "Νοσηλείες", "Νοσηλείες / 100.000 πληθυσμού", "Νοσηλείες / 100.000 πληθυσμού, με ηλικιακή στάθμιση",
  "Θάνατοι", "Θάνατοι / 100.000 πληθυσμού", "Θάνατοι / 100.000 πληθυσμού, με ηλικιακή στάθμιση")


# Extract the Vaccine Effectiveness estimates from the quasi-Poisson fits
getVE <- function(m) {
  res <- as.data.frame(100*(1 - exp(coef(summary(m))[paste0("vaccStatus",levels(dat$vaccStatus)[-1]),1:2] %*% rbind(1, c(0, 1.96, -1.96)))))
  colnames(res) <- c("est", "lo", "hi")
  rownames(res) <- paste0(1:3, "dose")
  res
}
VE <- lapply(m, getVE)


# Function to draw the Vaccine Effectiveness plot
out$VEplot <- function() { # width 8x6, pointsize 11
  par(family="Fira Sans", mar=c(5,7,2,11), oma=c(0,3,0,3))
  plot(0, type="n", bty="n", xlim=c(50,100), ylim=c(1,12), axes=F, ylab=NA, xlab="Αποτελεσματικότητα εμβολίου (%)", xaxs="i")
  with(VE$case, plotCI(est, 11:9, li=lo, ui=hi, err="x", add=TRUE, sfrac=0.00, lwd=2, pch=15, col=brewer.pal(4, "Spectral")[2:4]))
  with(VE$hosp, plotCI(est, 7:5, li=lo, ui=hi, err="x", add=TRUE, sfrac=0.00, lwd=2, pch=15, col=brewer.pal(4, "Spectral")[2:4]))
  with(VE$death, plotCI(est, 3:1, li=lo, ui=hi, err="x", add=TRUE, sfrac=0.00, lwd=2, pch=15, col=brewer.pal(4, "Spectral")[2:4]))
  axis(1)
  axis(2, at=c(11:9,7:5,3:1), labels=rep(c("1 δόση", "2 δόσεις", "3 δόσεις"),3), las=2, lwd=0)
  mtext("Έναντι εργαστηριακά επιβεβαιωμένης νόσου COVID-19", side=2, at=12, las=2, adj=0, line=2, font=2, cex=1.0)
  mtext("Έναντι νοσηλείας με COVID-19", side=2, at=8, las=2, adj=0, line=2, font=2, cex=1.0)
  mtext("Έναντι θανάτου με COVID-19", side=2, at=4, las=2, adj=0, line=2, font=2, cex=1.0)
  abline(v=100)
  axis(4, at=c(3:1,7:5,11:9), labels=with(do.call(rbind, VE), sprintf("%.1f (%.1f - %.1f)", est, lo, hi)), las=2, lwd=0)
  axis(4, at=12, labels="Αποτελεσματικότητα\nεμβολίου (%)", las=2, lwd=0)
}


# Function to draw all K-M plots
out$KMplot <- function() { # width 8x8, pointsize 11
  par(mfcol=c(3,2), family="Fira Sans", mar=c(5,5,4,1))
  plotAllKM("case", from="2021-2-1", to="2021-11-1", show=1:2, adj=FALSE, title = "Κρούσματα COVID-19, από 1/2/2021 έως 1/11/2021")
  plotAllKM("hosp", from="2021-2-1", to="2021-11-1", show=1:2, title = "Νοσηλείες COVID-19, από 1/2/2021 έως 1/11/2021")
  plotAllKM("death", from="2021-2-1", to="2021-11-1", show=1:2,  title = "Θάνατοι COVID-19, από 1/2/2021 έως 1/11/2021")
  plotAllKM("case", from="2021-11-2", adj=FALSE, title = "Κρούσματα COVID-19, από 1/11/2021 έως σήμερα")
  plotAllKM("hosp", from="2021-11-2", title = "Νοσηλείες COVID-19, από 1/11/2021 έως σήμερα")
  plotAllKM("death", from="2021-11-2", title = "Θάνατοι COVID-19, από 1/11/2021 έως σήμερα")
}



# Calculate various info "pieces" used in the .odt report
out$vaccDosesPerType <- sort(with(vacc, tapply(n, vaccType, sum))[vaccTypes[vaccTypes!="other"]], decreasing=TRUE)
out$saved <- lapply(predm, function(p) round(sum(p$x) + 1.96*c(0,-1,1)*sqrt(sum(p$var * p$x^2)) - sum(p$o)))
out$savedP <- list(
  death = do.call(sprintf, c("%s θανάτους (95%% CI: %s - %s)", as.list(out$saved$death))),
  hosp = do.call(sprintf, c("%s νοσηλείες (95%% CI: %s - %s)", as.list(out$saved$hosp))),
  case = do.call(sprintf, c("%s κρούσματα (95%% CI: %s - %s)", as.list(out$saved$case)))
)


# Set up output for odfWeave
library(odfWeave)

styleDefs <- getStyleDefs()
styleDefs$ArialCentered$fontName <- "Fira Sans"
styleDefs$ArialCentered$fontSize <- "10pt"
styleDefs$ArialCenteredBold$fontName <- "Fira Sans"
styleDefs$ArialCenteredBold$fontSize <- "10pt"
styleDefs$noBorder$verticalAlign <- "middle"
styleDefs$ArialRight <- styleDefs$ArialCentered
styleDefs$ArialRight$textAlign <- "right"
setStyleDefs(styleDefs)

# Style used in the rates tables
out$rtabstyle <- tableStyles(out$rtab1, header=colnames(out$rtab1))
out$rtabstyle$text[,1] <- "ArialRight"

# Pretty printer for numbers and percentages
out$pct1 <- function(x, y) sprintf("%.0f (%.1f%%)", x, x*100/y)



# Create the report, save everything and say goodbye
odfWeave("report-t.odt", "report.odt")


cat("Saving the output...\n")
save.image("intermediate.RData")

cat("Finished all analyses, have a nice day.\n")
