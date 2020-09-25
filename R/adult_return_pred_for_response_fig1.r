
## For generating predicted values for adult returns in Figure 1 of:
## "Faulkner, JR, BL Bellerud, DL Widener, SG Smith, and RW Zabel. 2020. Associations among fish length, dam passage history, and survival to adulthood in two at-risk species of Pacific salmon: response to comment."

## This script assumes there is a "data" directory containing all of the adult return data sets

# These packages need to be installed and loaded
library(lme4)


###-------------------------------
##           Data
###-------------------------------

### --- Read data files
cdat.ulgr <- read.csv("data/chinook_ULGR_SAR_data.csv", stringsAsFactors=F)
cdat.lgr <- read.csv("data/chinook_LGR_SAR_data.csv", stringsAsFactors=F)
sdat.ulgr <- read.csv("data/steelhead_ULGR_SAR_data.csv", stringsAsFactors=F)
sdat.lgr <- read.csv("data/steelhead_LGR_SAR_data.csv", stringsAsFactors=F)


## --- Make new variables
## make year a factor
cdat.ulgr$fyear <- as.factor(cdat.ulgr$year)
cdat.lgr$fyear <- as.factor(cdat.lgr$year)
sdat.ulgr$fyear <- as.factor(sdat.ulgr$year)
sdat.lgr$fyear <- as.factor(sdat.lgr$year)



#### ------------------------------------------------
##      Generate Predictions for Chinook
##      Tagged Upstream of Lower Granite
#### ------------------------------------------------

# Refit best model from original analysis using categorical number of bypasses

## --- Make categorical bypass variable with 5+ category
cdat.ulgr$cnbypass <- as.character(cdat.ulgr$nbypass)
cdat.ulgr$cat.nbypass.5plus <- as.character(cdat.ulgr$nbypass)
cdat.ulgr$cat.nbypass.5plus[cdat.ulgr$nbypass>=5] <- "5"

## --- Fit model
ch.fit.ulgr <- glmer(adult_return ~ zlength + cat.nbypass.5plus + clipped + rel_site2 + zlast_day + last_site + (1 | fyear), data=cdat.ulgr, family=binomial)
summary(ch.fit.ulgr)
## Refit for convergence
tmp.refit <- getME(ch.fit.ulgr,c("theta","fixef"))
ch.fit.ulgr2 <- update(ch.fit.ulgr,start=tmp.refit,control=glmerControl(optCtrl=list(maxfun=2e5)) )
summary(ch.fit.ulgr2)

## -- Make prediction data set
## get number of fish last seen at trawl
ch.lstab.ulgr = table(cdat.ulgr$last_site)
ch.prop.twx.ulgr = ch.lstab.ulgr[2]/sum(ch.lstab.ulgr)
## create design matrix
# columns: ones, zlength=0, cat.nbypass.5plus.1, cat.nbypass.5plus.2 ,cat.nbypass.5plus.3,cat.nbypass.5plus.4,
#    cat.nbypass.5plus.5p, clipped=no, rel_site2=CLSNK, zlast_day=0, last_site
# note last site is weighted average between BON and TWX based on number of fish, and year effect is zero (mean) 
ch.xmat.ulgr = cbind(rep(1,6), 0, c(0,1,0,0,0,0), c(0,0,1,0,0,0),c(0,0,0,1,0,0),c(0,0,0,0,1,0), c(0,0,0,0,0,1), rep(1,6), rep(0,6), rep(0,6), rep(ch.prop.twx.ulgr,6))  

## -- Extract fixed effect parameters
ch.betav.ulgr = summary(ch.fit.ulgr2)$coef[,1]

## -- Calculate predicted values and CI's
ch.pred.ulgr <- plogis(ch.xmat.ulgr%*%ch.betav.ulgr)   
ch.vc.ulgr <- ch.xmat.ulgr%*%vcov(ch.fit.ulgr2)%*%t(ch.xmat.ulgr)  #variance of fixed effects
ch.cil.ulgr = plogis(qlogis(ch.pred.ulgr) - 1.96*diag(sqrt(ch.vc.ulgr)))
ch.ciu.ulgr = plogis(qlogis(ch.pred.ulgr) + 1.96*diag(sqrt(ch.vc.ulgr)))

## -- Make data frame with outputs
ch.tab5.ulgr = table(cdat.ulgr$cat.nbypass.5plus, cdat.ulgr$adult_return)
ch.ulgr.out <- data.frame(nbypass=0:5, nfish=ch.tab5.ulgr[,1]+ch.tab5.ulgr[,2], return=ch.tab5.ulgr[,2], noreturn=ch.tab5.ulgr[,1],pred=ch.pred.ulgr, cil=ch.cil.ulgr, ciu=ch.ciu.ulgr)




#### ------------------------------------------------
##      Generate Predictions for  
##      Chinook Tagged at Lower Granite
#### ------------------------------------------------

## --- Make categorical bypass variable with 5+ category
## first add 1 to nbypass (consider all lgr fish as bypassed)
cdat.lgr$nbypass1 = cdat.lgr$nbypass + 1
cdat.lgr$cat.nbypass.5plus <- as.character(cdat.lgr$nbypass1)
cdat.lgr$cat.nbypass.5plus[cdat.lgr$nbypass1>=5] <- "5"

## --- Fit model
ch.fit.lgr <- glmer(adult_return ~ zlength + cat.nbypass.5plus + clipped + zlast_day + (zlast_day | fyear), data=cdat.lgr, family=binomial)
summary(ch.fit.lgr)


## -- Make prediction data set
## create design matrix
# columns: ones, zlength=0, cat.nbypass.5plus.2 ,cat.nbypass.5plus.3,cat.nbypass.5plus.4,
#    cat.nbypass.5plus.5p, clipped=no, zlast_day=0
# note year effects are zero (mean) 
ch.xmat.lgr = cbind(rep(1,5), 0, c(0,1,0,0,0), c(0,0,1,0,0),c(0,0,0,1,0),c(0,0,0,0,1),  rep(1,5), rep(0,5))

## -- Extract fixed effect parameters
ch.betav.lgr = summary(ch.fit.lgr)$coef[,1]

## -- Calculate predicted values and CI's
ch.pred.lgr <- plogis(ch.xmat.lgr%*%ch.betav.lgr)
ch.vc.lgr <- ch.xmat.lgr%*%vcov(ch.fit.lgr)%*%t(ch.xmat.lgr)
ch.cil.lgr = plogis(qlogis(ch.pred.lgr) - 1.96*diag(sqrt(ch.vc.lgr)))
ch.ciu.lgr = plogis(qlogis(ch.pred.lgr) + 1.96*diag(sqrt(ch.vc.lgr)))

## -- Make data frame with outputs
ch.tab5.lgr = table(cdat.lgr$cat.nbypass.5plus, cdat.lgr$adult_return)
ch.lgr.out <- data.frame(nbypass=1:5,  nfish=ch.tab5.lgr[,1]+ch.tab5.lgr[,2], return=ch.tab5.lgr[,2], noreturn=ch.tab5.lgr[,1],pred=ch.pred.lgr, cil=ch.cil.lgr, ciu=ch.ciu.lgr)






#### ------------------------------------------------
##      Generate Predictions for Steelhead
##      Tagged Upstream of Lower Granite
#### ------------------------------------------------

# Refit best model from original analysis using categorical number of bypasses

## --- Make categorical bypass variable with 5+ category
sdat.ulgr$cnbypass <- as.character(sdat.ulgr$nbypass)
sdat.ulgr$cat.nbypass.5plus <- as.character(sdat.ulgr$nbypass)
sdat.ulgr$cat.nbypass.5plus[sdat.ulgr$nbypass>=5] <- "5"

## --- Fit model
st.fit.ulgr <- glmer(adult_return ~ zlength + cat.nbypass.5plus + clipped + rel_site2 + zlast_day + zlast_day_sq + last_site + (1 | fyear), data=sdat.ulgr, family=binomial)
summary(st.fit.ulgr)

## -- Make prediction data set
## get number of fish last seen at trawl
st.lstab.ulgr = table(sdat.ulgr$last_site)
st.prop.twx.ulgr = st.lstab.ulgr[2]/sum(st.lstab.ulgr)
## create design matrix
# columns: ones, zlength=0, cat.nbypass.5plus.1, cat.nbypass.5plus.2 ,cat.nbypass.5plus.3,cat.nbypass.5plus.4,
#    cat.nbypass.5plus.5p, clipped=no, rel_site2=CLSNK, zlast_day=0,zlast_day_sq=0, last_site
# note last site is weighted average between BON and TWX based on number of fish, and year effect is zero (mean) 
st.xmat.ulgr = cbind(rep(1,6), 0, c(0,1,0,0,0,0), c(0,0,1,0,0,0),c(0,0,0,1,0,0),c(0,0,0,0,1,0), c(0,0,0,0,0,1), rep(1,6), rep(0,6), rep(0,6), rep(0,6), rep(st.prop.twx.ulgr,6))  

## -- Extract fixed effect parameters
st.betav.ulgr = summary(st.fit.ulgr)$coef[,1]

## -- Calculate predicted values and CI's
st.pred.ulgr <- plogis(st.xmat.ulgr%*%st.betav.ulgr)   
st.vc.ulgr <- st.xmat.ulgr%*%vcov(st.fit.ulgr)%*%t(st.xmat.ulgr)  #variance of fixed effects
st.cil.ulgr = plogis(qlogis(st.pred.ulgr) - 1.96*diag(sqrt(st.vc.ulgr)))
st.ciu.ulgr = plogis(qlogis(st.pred.ulgr) + 1.96*diag(sqrt(st.vc.ulgr)))

## -- Make data frame with outputs
st.tab5.ulgr = table(sdat.ulgr$cat.nbypass.5plus, sdat.ulgr$adult_return)
st.ulgr.out <- data.frame(nbypass=0:5, nfish=st.tab5.ulgr[,1]+st.tab5.ulgr[,2], return=st.tab5.ulgr[,2], noreturn=st.tab5.ulgr[,1],pred=st.pred.ulgr, cil=st.cil.ulgr, ciu=st.ciu.ulgr)




#### ------------------------------------------------
##      Generate Predictions for  
##      Steelhead Tagged at Lower Granite
#### ------------------------------------------------

## --- Make categorical bypass variable with 5+ category
## first add 1 to nbypass (consider all lgr fish as bypassed)
sdat.lgr$nbypass1 = sdat.lgr$nbypass + 1
sdat.lgr$cat.nbypass.5plus <- as.character(sdat.lgr$nbypass1)
sdat.lgr$cat.nbypass.5plus[sdat.lgr$nbypass1>=5] <- "5"

## --- Fit model
st.fit.lgr <- glmer(adult_return ~ zlength + cat.nbypass.5plus + clipped + zlast_day + (zlast_day | fyear), data=sdat.lgr, family=binomial)
summary(st.fit.lgr)


## -- Make prediction data set
## create design matrix
# columns: ones, zlength=0, cat.nbypass.5plus.2 ,cat.nbypass.5plus.3,cat.nbypass.5plus.4,
#    cat.nbypass.5plus.5p, clipped=no, zlast_day=0
# note year effects are zero (mean) 
st.xmat.lgr = cbind(rep(1,5), 0, c(0,1,0,0,0), c(0,0,1,0,0),c(0,0,0,1,0),c(0,0,0,0,1),  rep(1,5), rep(0,5))

## -- Extract fixed effect parameters
st.betav.lgr = summary(st.fit.lgr)$coef[,1]

## -- Calculate predicted values and CI's
st.pred.lgr <- plogis(st.xmat.lgr%*%st.betav.lgr)
st.vc.lgr <- st.xmat.lgr%*%vcov(st.fit.lgr)%*%t(st.xmat.lgr)
st.cil.lgr = plogis(qlogis(st.pred.lgr) - 1.96*diag(sqrt(st.vc.lgr)))
st.ciu.lgr = plogis(qlogis(st.pred.lgr) + 1.96*diag(sqrt(st.vc.lgr)))

## -- Make data frame with outputs
st.tab5.lgr = table(sdat.lgr$cat.nbypass.5plus, sdat.lgr$adult_return)
st.lgr.out <- data.frame(nbypass=1:5,  nfish=st.tab5.lgr[,1]+st.tab5.lgr[,2], return=st.tab5.lgr[,2], noreturn=st.tab5.lgr[,1],pred=st.pred.lgr, cil=st.cil.lgr, ciu=st.ciu.lgr)






### --------------------------------------------------
##          Make Plot
### --------------------------------------------------

poff = 0.12
colv <- c("gray10", "gray50", "gray30")
cylim = c(0,.055) #c(0,.04)
sylim = c(0,.11)
byplab = c("0","1","2","3","4","5-7")

cyat = c(0.00, 0.01, 0.02, 0.03, 0.04, 0.05)
cylabs = c("0.00", "", "0.02", "", "0.04", "")
syat = c(0.00,0.02, 0.04,0.06, 0.08,0.10, 0.12) #
sylabs = c("0.00", "", "0.04", "", "0.08", "", "0.12")

pdf("fig1_return_probabilities.pdf", height=5, width=9)
 par(mfrow=c(1, 2), oma=c(0.25,1,0,0), mar=c(4,3,2.5,1), font.axis=2 )

 plot(c(0:4,6), ch.ulgr.out$pred, ylim=cylim, xlim=c(-0.5, 6.5), type="n", xaxt="n", yaxt="n", xlab="", ylab="", main="Chinook Salmon")
 axis(side=1, at=c(0:4,6), labels=byplab)
 axis(side=2, at=cyat, labels=cylabs)
 abline(h=ch.ulgr.out$pred[1], lty=5, col=colv[3])
 abline(h=c(ch.ulgr.out$cil[1], ch.ulgr.out$ciu[1]), lty=2, col=colv[3])
 points(c(0:4,6)-poff, ch.ulgr.out$pred, pch=16, cex=1.4, col=colv[1])
 arrows(x0=c(0:4,6)-poff, y0=ch.ulgr.out$cil, x1=c(0:4,6)-poff, y1=ch.ulgr.out$ciu, angle=90, code=3, length=.05, col=colv[1], lwd=2)
 points(c(1:4,6)+poff, ch.lgr.out$pred, pch=16, cex=1.4, col=colv[2])
 arrows(x0=c(1:4,6)+poff, y0=ch.lgr.out$cil, x1=c(1:4,6)+poff, y1=ch.lgr.out$ciu, angle=90, code=3, length=.05, col=colv[2], lwd=2)
 legend("topleft", legend=c("Tagged ULGR", "Tagged LGR"), pch=c(16,16), lty=c(1,1), col=colv[1:2], bty="n", cex=1, pt.cex=1.4, lwd=2)

 plot(c(0:4,6), st.ulgr.out$pred, ylim=sylim, xlim=c(-0.5, 6.5), type="n", xaxt="n", yaxt="n", xlab="", ylab="", main="Steelhead")
 axis(side=1, at=c(0:4,6), labels=byplab)
 axis(side=2, at=syat, labels=sylabs)

 abline(h=st.ulgr.out$pred[1], lty=5, col=colv[3])
 abline(h=c(st.ulgr.out$cil[1], st.ulgr.out$ciu[1]), lty=2, col=colv[3])
 points(c(0:4,6)-poff, st.ulgr.out$pred, pch=16, cex=1.4, col=colv[1])
 arrows(x0=c(0:4,6)-poff, y0=st.ulgr.out$cil, x1=c(0:4,6)-poff, y1=st.ulgr.out$ciu, angle=90, code=3, length=.05, col=colv[1], lwd=2)
 points(c(1:4,6)+poff, st.lgr.out$pred, pch=16, cex=1.4, col=colv[2])
 arrows(x0=c(1:4,6)+poff, y0=st.lgr.out$cil, x1=c(1:4,6)+poff, y1=st.lgr.out$ciu, angle=90, code=3, length=.05, col=colv[2], lwd=2)

 mtext(side=1, line=-0.75, outer=T, text="Number of bypass events", font=2)
 mtext(side=2, line=0, outer=T, text="Probability of adult return", font=2)

dev.off()










