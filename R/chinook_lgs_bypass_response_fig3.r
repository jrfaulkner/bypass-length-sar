

## R code for generating predicted bypass probabilities used in Figure 3 of:
### "Faulkner, J.R., B.L. Bellerud, D.L. Widener, S.G. Smith, and R.W. Zabel.
###   Associations among fish length, dam passage history, and survival to adulthood
###   in two at-risk species of Pacific salmon: response to comment." 

## This script assumes there is a "data" directory with the chinook bypass data


library(lme4)



###-------------------------------
##         Data Processing
###-------------------------------

### --- Read data files
cdat <- read.csv("data/chinook_bypass_data.csv", stringsAsFactors=F)


## --- Make new variables
# make year a factor
cdat$fyear <- as.character(cdat$year)

## --- Get quantiles of length and day
# quantiles for standardized length
qvals = quantile(cdat$zlength, probs=c(.01, .05, .95, .99))

# quantiles for day of passage at LGR
yrset = 2000:2014
tmp.q1 = aggregate(cdat$day_lgr, by=list(year=cdat$fyear), quantile,probs=c(.01))$x
tmp.q2 = aggregate(cdat$day_lgr, by=list(year=cdat$fyear), quantile,probs=c(.99))$x
dayqmat = data.frame(year=yrset, q01=tmp.q1, q99=tmp.q2) 



### -----------------------------------------------------
##      Model Fitting and Parameter
### -----------------------------------------------------

## Fit best model for Little Goose Dam from original paper
modfit = glmer(bypass_LGS ~ zlength + clipped + rel_site2 + zday_lgr + zday_lgr_sq + (zday_lgr + zday_lgr_sq | fyear), family=binomial, data=cdat)


gsdat = coef(modfit)$fyear
names(gsdat)[c(1,3)] = c("intercept", "clippedb.no")



### ------------------------------------
##      Plot of Subset of Years
### ------------------------------------


syrset = c(2004, 2005, 2007, 2010, 2011, 2014)


xrngz = c(-2,3)
xrngb = mean(cdat$day_lgr) + sd(cdat$day_lgr)*xrngz


sdlabs = c(110, 120, 130, 140, 150) #round(sdtix*10.67859 + 126.01) 
sdtix = (c(110, 120, 130, 140, 150)-126.01)/10.6786  #c(-2, -1, 0, 1, 2, 3)

pdf("Fig_3_chinook_bypass_predictions.pdf", height=7, width=7 )
par(mfrow=c(3,2), mar=c(2,3.5,0.5,1), oma=c(2,1,0,0), cex.axis=1.2, font.axis=2 )
for (i in 1:length(syrset)) {
  tmp.yrid = which(dayqmat$year==syrset[i])
  tmp.dseq = seq(max(xrngb[1], dayqmat[tmp.yrid,2]), min(xrngb[2], dayqmat[tmp.yrid,3]), length=200)
  tmp.zdseq = (tmp.dseq-mean(cdat$day_lgr))/sd(cdat$day_lgr)
  tmp.zdseq2 = tmp.zdseq^2

  tmp.pred = plogis(gsdat$intercept[tmp.yrid] + gsdat$clippedb.no[tmp.yrid] + gsdat$zday_lgr[tmp.yrid]*tmp.zdseq + gsdat$zday_lgr_sq[tmp.yrid]*tmp.zdseq2)
  tmp.pred.q1 = plogis(gsdat$intercept[tmp.yrid] + gsdat$zlength[tmp.yrid]*qvals[1] + gsdat$clippedb.no[tmp.yrid] + gsdat$zday_lgr[tmp.yrid]*tmp.zdseq + gsdat$zday_lgr_sq[tmp.yrid]*tmp.zdseq2)
    tmp.pred.q4 = plogis(gsdat$intercept[tmp.yrid] + gsdat$zlength[tmp.yrid]*qvals[4] + gsdat$clippedb.no[tmp.yrid] + gsdat$zday_lgr[tmp.yrid]*tmp.zdseq + gsdat$zday_lgr_sq[tmp.yrid]*tmp.zdseq2)

  plot(tmp.zdseq, tmp.pred, type="l", xlab="", ylab="", main="", ylim=c(0,1), xlim=xrngz, lwd=2, xaxt="n" )
  if (i %in% 1:4) axis(side=1, at=sdtix, labels=rep("",5) )
  if (i %in% 5:6) axis(side=1, at=sdtix, labels=sdlabs )
  lines(tmp.zdseq, tmp.pred.q1, lty=2, lwd=2)
  lines(tmp.zdseq, tmp.pred.q4, lty=4, lwd=2)
  
  text(x=-1.75, y=.97, syrset[i], font=2, cex=1.5)
  #if (i==2) legend("bottomleft", title="Length", legend=c("Mean", "1%", "99%"), lty=c(1,2,4), lwd=2, bty="n", cex=1.2)
  if (i==2) legend("bottomleft", legend=c("Mean length", "Small fish (1%)", "Large fish (99%)"), lty=c(1,2,4), lwd=2, bty="n", cex=1.2)
  
}
  mtext(side=1, line=0.9, outer=T, text="Day of year", font=2)
  mtext(side=2, line=-0.2, outer=T, text="Bypass probability", font=2)

dev.off()








