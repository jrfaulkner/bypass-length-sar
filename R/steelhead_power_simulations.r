
### R Code for calculating power to detect bypass effects for steelhead at
### various effect sizes and sample sizes.
### This analysis is part of 
### "Faulkner, J.R., B.L. Bellerud, D.L. Widener, S.G. Smith, and R.W. Zabel.
###   Associations among fish length, dam passage history, and survival to adulthood
###   in two at-risk species of Pacific salmon: response to comment." 

## This script assumes there is a "data" directory with the SAR data and 
## an "output" directory for writing output files

library(lme4)


### --------------------------------
##          Data
### --------------------------------

## -- Tagged upstream of LGR --------
# data used for SAR model fitting
fsdat.ulgr = read.csv("data/steelhead_ULGR_SAR_data.csv", stringsAsFactors=F) 
fsdat.ulgr$fyear = factor(fsdat.ulgr$year)

yrset = 2004:2014
nyear = length(yrset)

# Table for number of bypass events by year
byptab.ulgr = table(fsdat.ulgr$year, fsdat.ulgr$nbypass)
# Create data frame from table
sdat1.ulgr = as.data.frame(byptab.ulgr)
names(sdat1.ulgr) = c("year", "nbypass", "nfish")
sdat1.ulgr = sdat1.ulgr[order(sdat1.ulgr$year, sdat1.ulgr$nbypass), ]


## -- Tagged at LGR -----------
# data used for SAR model fitting
fsdat.lgr = read.csv("data/steelhead_LGR_SAR_data.csv", stringsAsFactors=F) 
fsdat.lgr$fyear = factor(fsdat.lgr$year)


# Table for number of bypass events by year
byptab.lgr = table(fsdat.lgr$year, fsdat.lgr$nbypass)
# Create data frame from table
sdat1.lgr = as.data.frame(byptab.lgr)
names(sdat1.lgr) = c("year", "nbypass", "nfish")
sdat1.lgr = sdat1.lgr[order(sdat1.lgr$year, sdat1.lgr$nbypass), ]


## ----------------------------------------------------
##      Initial Model Fitting for Parameter Values 
## ----------------------------------------------------

## ULGR - model for intercept and year effects
tfit.ulgr = glmer(adult_return ~ nbypass + (1 | fyear), data=fsdat.ulgr, family=binomial)
# extract random year effects
ranef.yr.ulgr = as.vector(as.matrix(ranef(tfit.ulgr)$fyear))

## LGR - model for intercept and year effects
tfit.lgr = glmer(adult_return ~ nbypass + (1 | fyear), data=fsdat.lgr, family=binomial)
# extract random year effects
ranef.yr.lgr = as.vector(as.matrix(ranef(tfit.lgr)$fyear))



## ----------------------------------------------------
##      Power Simulation for 
##      Fish Tagged Upstream of Lower Granite Dam
## ----------------------------------------------------

## Set parameters
tmp.int.ulgr = as.numeric(fixef(tfit.ulgr)[1])  #intercept
tmp.yreff.ulgr = rep(ranef.yr.ulgr, rep(8, nyear)) #year effects 
tmp.nbyp.alt1 = -0.10   # alternative value 1
tmp.nbyp.alt2 = -0.15   # alternative value 2
tmp.nbyp.css1 = -0.127  # 2010 value
tmp.nbyp.css2 = -0.116  # 2016 value

## Calculate return probabilities
sdat1.ulgr$nu.alt1 = tmp.int.ulgr + tmp.yreff.ulgr + tmp.nbyp.alt1*sdat1.ulgr$nbypass
sdat1.ulgr$psar.alt1 = plogis(sdat1.ulgr$nu.alt1)
sdat1.ulgr$nu.alt2 = tmp.int.ulgr + tmp.yreff.ulgr + tmp.nbyp.alt2*sdat1.ulgr$nbypass
sdat1.ulgr$psar.alt2 = plogis(sdat1.ulgr$nu.alt2)
sdat1.ulgr$nu.css1 = tmp.int.ulgr + tmp.yreff.ulgr + tmp.nbyp.css1*sdat1.ulgr$nbypass
sdat1.ulgr$psar.css1 = plogis(sdat1.ulgr$nu.css1)
sdat1.ulgr$nu.css2 = tmp.int.ulgr + tmp.yreff.ulgr + tmp.nbyp.css2*sdat1.ulgr$nbypass
sdat1.ulgr$psar.css2 = plogis(sdat1.ulgr$nu.css2)
## drop groupings with zero fish
sdat.ulgr = sdat1.ulgr[sdat1.ulgr$nfish!=0, ]
sdat.ulgr$fyear = as.character(sdat.ulgr$year)

nobs.ulgr = nrow(sdat.ulgr) #number of observation groupings

nsim = 1000  #number of simulations
## set up matrices to collect outputs
simout.ulgr.css1 = matrix(0, nsim, 5)  #parm, se, pval, aic.yronly, aic.full
simout.ulgr.css2 = matrix(0, nsim, 5)
simout.ulgr.alt1 = matrix(0, nsim, 5)
simout.ulgr.alt2 = matrix(0, nsim, 5)

## Run simulations
for (k in 1:nsim) {
  cat("running sim", k, "of", nsim, "\n")

  ## Bypass effect = Tuomikoski et al 2010 value 
  # generate random number of returning fish
  set.seed(k*10)
  tmp.return.css1 = rbinom(nobs.ulgr, size=sdat.ulgr$nfish, prob=sdat.ulgr$psar.css1)
  tmp.noreturn.css1 = as.integer(sdat.ulgr$nfish - tmp.return.css1)
  tmp.adult_return.css1 = cbind(tmp.return.css1, tmp.noreturn.css1)
  # fit models with and without nbypass
  tmp.fit.css1 = glmer(tmp.adult_return.css1 ~ nbypass + (1 | fyear), data=sdat.ulgr, family=binomial)
  tmp.fit.css1.yr = glmer(tmp.adult_return.css1 ~ (1 | fyear), data=sdat.ulgr, family=binomial)
  tmp.sum.css1 = summary(tmp.fit.css1)
  # collect outputs 
  simout.ulgr.css1[k,] = c(as.vector(tmp.sum.css1$coef[2, -3]), AIC(tmp.fit.css1.yr), AIC(tmp.fit.css1)) 

  ## Bypass effect = McCann et al 2016 value 
  # generate random number of returning fish
  set.seed(k*10)
  tmp.return.css2 = rbinom(nobs.ulgr, size=sdat.ulgr$nfish, prob=sdat.ulgr$psar.css2)
  tmp.noreturn.css2 = as.integer(sdat.ulgr$nfish - tmp.return.css2)
  tmp.adult_return.css2 = cbind(tmp.return.css2, tmp.noreturn.css2)
  # fit models with and without nbypass
  tmp.fit.css2 = glmer(tmp.adult_return.css2 ~ nbypass + (1 | fyear), data=sdat.ulgr, family=binomial)
  tmp.fit.css2.yr = glmer(tmp.adult_return.css2 ~ (1 | fyear), data=sdat.ulgr, family=binomial)
  tmp.sum.css2 = summary(tmp.fit.css2)
  # collect outputs 
  simout.ulgr.css2[k,] = c(as.vector(tmp.sum.css2$coef[2, -3]), AIC(tmp.fit.css2.yr), AIC(tmp.fit.css2))

  ## Bypass effect = Alternative value 1 
  # generate random number of returning fish
  set.seed(k*10)
  tmp.return.alt1 = rbinom(nobs.ulgr, size=sdat.ulgr$nfish, prob=sdat.ulgr$psar.alt1)
  tmp.noreturn.alt1 = as.integer(sdat.ulgr$nfish - tmp.return.alt1)
  tmp.adult_return.alt1 = cbind(tmp.return.alt1, tmp.noreturn.alt1)
  # fit models with and without nbypass
  tmp.fit.alt1 = glmer(tmp.adult_return.alt1 ~ nbypass + (1 | fyear), data=sdat.ulgr, family=binomial)
  tmp.fit.alt1.yr = glmer(tmp.adult_return.alt1 ~ (1 | fyear), data=sdat.ulgr, family=binomial)
  tmp.sum.alt1 = summary(tmp.fit.alt1)
  # collect outputs 
  simout.ulgr.alt1[k,] = c(as.vector(tmp.sum.alt1$coef[2, -3]), AIC(tmp.fit.alt1.yr), AIC(tmp.fit.alt1))

  ## Bypass effect = Alternative value 2 
  # generate random number of returning fish
  set.seed(k*10)
  tmp.return.alt2 = rbinom(nobs.ulgr, size=sdat.ulgr$nfish, prob=sdat.ulgr$psar.alt2)
  tmp.noreturn.alt2 = as.integer(sdat.ulgr$nfish - tmp.return.alt2)
  tmp.adult_return.alt2 = cbind(tmp.return.alt2, tmp.noreturn.alt2)
  # fit models with and without nbypass
  tmp.fit.alt2 = glmer(tmp.adult_return.alt2 ~ nbypass + (1 | fyear), data=sdat.ulgr, family=binomial)
  tmp.fit.alt2.yr = glmer(tmp.adult_return.alt2 ~  (1 | fyear), data=sdat.ulgr, family=binomial)
  tmp.sum.alt2 = summary(tmp.fit.alt2)
  # collect outputs 
  simout.ulgr.alt2[k,] = c(as.vector(tmp.sum.alt2$coef[2, -3]), AIC(tmp.fit.alt2.yr), AIC(tmp.fit.alt2))

}

## --- Calculate Power for Various Criteria ----
## -- alpha = 0.05
power.ulgr.css1 = sum(simout.ulgr.css1[,3]<=0.05)/nsim
power.ulgr.css2 = sum(simout.ulgr.css2[,3]<=0.05)/nsim
power.ulgr.alt1 = sum(simout.ulgr.alt1[,3]<=0.05)/nsim
power.ulgr.alt2 = sum(simout.ulgr.alt2[,3]<=0.05)/nsim

## -- alpha = 0.10
power.ulgr.css1.10 = sum(simout.ulgr.css1[,3]<=0.1)/nsim
power.ulgr.css2.10 = sum(simout.ulgr.css2[,3]<=0.1)/nsim
power.ulgr.alt1.10 = sum(simout.ulgr.alt1[,3]<=0.1)/nsim
power.ulgr.alt2.10 = sum(simout.ulgr.alt2[,3]<=0.1)/nsim

## -- delta AIC >= 2
power.ulgr.css1.A2 = sum((simout.ulgr.css1[,4]-simout.ulgr.css1[,5])>=2)/nsim
power.ulgr.css2.A2 = sum((simout.ulgr.css2[,4]-simout.ulgr.css2[,5])>=2)/nsim
power.ulgr.alt1.A2 = sum((simout.ulgr.alt1[,4]-simout.ulgr.alt1[,5])>=2)/nsim
power.ulgr.alt2.A2 = sum((simout.ulgr.alt2[,4]-simout.ulgr.alt2[,5])>=2)/nsim

## -- delta AIC >= 0
power.ulgr.css1.A0 = sum((simout.ulgr.css1[,4]-simout.ulgr.css1[,5])>=0)/nsim
power.ulgr.css2.A0 = sum((simout.ulgr.css2[,4]-simout.ulgr.css2[,5])>=0)/nsim
power.ulgr.alt1.A0 = sum((simout.ulgr.alt1[,4]-simout.ulgr.alt1[,5])>=0)/nsim
power.ulgr.alt2.A0 = sum((simout.ulgr.alt2[,4]-simout.ulgr.alt2[,5])>=0)/nsim


## Table for results
simsum.ulgr = data.frame(species="Steelhead", tagsite="ULGR", matrix(0,4,5))
names(simsum.ulgr)[3:7] = c("parameter", "alpha_05", "alpha_10", "AIC2", "AIC0")
simsum.ulgr$parameter = c(tmp.nbyp.css1,tmp.nbyp.css2, tmp.nbyp.alt1,tmp.nbyp.alt2)
simsum.ulgr$alpha_05 = c(power.ulgr.css1, power.ulgr.css2, power.ulgr.alt1, power.ulgr.alt2)
simsum.ulgr$alpha_10 = c(power.ulgr.css1.10, power.ulgr.css2.10, power.ulgr.alt1.10, power.ulgr.alt2.10)
simsum.ulgr$AIC2 = c(power.ulgr.css1.A2, power.ulgr.css2.A2, power.ulgr.alt1.A2, power.ulgr.alt2.A2)
simsum.ulgr$AIC0 = c(power.ulgr.css1.A0, power.ulgr.css2.A0, power.ulgr.alt1.A0, power.ulgr.alt2.A0)
simsum.ulgr = simsum.ulgr[order(simsum.ulgr$parameter, decreasing=T), ]



 

## --- Try with fixed effects for year --- 
# NOTE: results for this were referenced in the response to comment but not directly reported
# This was tested as an explanation for lower power reported by Storch et al
nsim = 1000
simout.ulgr.css1.fixef = matrix(0, nsim, 3)

for (k in 1:nsim) {
  cat("running sim", k, "of", nsim, "\n")
  
  set.seed(k*10)
  tmp.return.css = rbinom(nobs.ulgr, size=sdat.ulgr$nfish, prob=sdat.ulgr$psar.css1)
  tmp.noreturn.css = as.integer(sdat.ulgr$nfish - tmp.return.css)
  tmp.adult_return.css = cbind(tmp.return.css, tmp.noreturn.css)

  tmp.fit.css = glm(tmp.adult_return.css ~ nbypass + fyear, data=sdat.ulgr, family=binomial)
  tmp.sum.css = summary(tmp.fit.css)
  simout.ulgr.css1.fixef[k,] = as.vector(tmp.sum.css$coef[2, -3])

}

power.ulgr.css1.fixef = sum(simout.ulgr.css1.fixef[,3]<=0.05)/nsim
power.ulgr.css1.fixef





## --------------------------------------------
##      Power Simulation for 
##      Fish Tagged at Lower Granite Dam
## --------------------------------------------

## Set parameters
tmp.int.lgr = as.numeric(fixef(tfit.lgr)[1])  #intercept
tmp.yreff.lgr = rep(ranef.yr.lgr, rep(7, nyear)) #year effects 
tmp.nbyp.alt1 = -0.10   # alternative value 1
tmp.nbyp.alt2 = -0.15   # alternative value 2
tmp.nbyp.css1 = -0.127  # 2010 value
tmp.nbyp.css2 = -0.116  # 2016 value

## Calculate return probabilities
sdat1.lgr$nu.alt1 = tmp.int.lgr + tmp.yreff.lgr + tmp.nbyp.alt1*sdat1.lgr$nbypass
sdat1.lgr$psar.alt1 = plogis(sdat1.lgr$nu.alt1)
sdat1.lgr$nu.alt2 = tmp.int.lgr + tmp.yreff.lgr + tmp.nbyp.alt2*sdat1.lgr$nbypass
sdat1.lgr$psar.alt2 = plogis(sdat1.lgr$nu.alt2)
sdat1.lgr$nu.css1 = tmp.int.lgr + tmp.yreff.lgr + tmp.nbyp.css1*sdat1.lgr$nbypass
sdat1.lgr$psar.css1 = plogis(sdat1.lgr$nu.css1)
sdat1.lgr$nu.css2 = tmp.int.lgr + tmp.yreff.lgr + tmp.nbyp.css2*sdat1.lgr$nbypass
sdat1.lgr$psar.css2 = plogis(sdat1.lgr$nu.css2)

## drop groupings with zero fish
sdat.lgr = sdat1.lgr[sdat1.lgr$nfish!=0, ]
sdat.lgr$fyear = as.character(sdat.lgr$year)

nobs.lgr = nrow(sdat.lgr)  #number of observation groupings

nsim = 1000    #number of simulations
## set up matrices to collect outputs
simout.lgr.css1 = matrix(0, nsim, 5)  #parm, se, pval, aic.yronly, aic.full
simout.lgr.css2 = matrix(0, nsim, 5)
simout.lgr.alt1 = matrix(0, nsim, 5)
simout.lgr.alt2 = matrix(0, nsim, 5)

## Run simulations
for (k in 1:nsim) {
  cat("running sim", k, "of", nsim, "\n")

  ## Bypass effect = Tuomikoski et al 2010 value 
  # generate random number of returning fish
  set.seed(k*10)
  tmp.return.css1 = rbinom(nobs.lgr, size=sdat.lgr$nfish, prob=sdat.lgr$psar.css1)
  tmp.noreturn.css1 = as.integer(sdat.lgr$nfish - tmp.return.css1)
  tmp.adult_return.css1 = cbind(tmp.return.css1, tmp.noreturn.css1)
  # fit models with and without nbypass
  tmp.fit.css1 = glmer(tmp.adult_return.css1 ~ nbypass + (1 | fyear), data=sdat.lgr, family=binomial)
  tmp.fit.css1.yr = glmer(tmp.adult_return.css1 ~ (1 | fyear), data=sdat.lgr, family=binomial)
  tmp.sum.css1 = summary(tmp.fit.css1)
  # collect outputs
  simout.lgr.css1[k,] = c(as.vector(tmp.sum.css1$coef[2, -3]), AIC(tmp.fit.css1.yr), AIC(tmp.fit.css1)) 

  ## Bypass effect = McCann et al 2016 value 
  # generate random number of returning fish
  set.seed(k*10)
  tmp.return.css2 = rbinom(nobs.lgr, size=sdat.lgr$nfish, prob=sdat.lgr$psar.css2)
  tmp.noreturn.css2 = as.integer(sdat.lgr$nfish - tmp.return.css2)
  tmp.adult_return.css2 = cbind(tmp.return.css2, tmp.noreturn.css2)
  # fit models with and without nbypass
  tmp.fit.css2 = glmer(tmp.adult_return.css2 ~ nbypass + (1 | fyear), data=sdat.lgr, family=binomial)
  tmp.fit.css2.yr = glmer(tmp.adult_return.css2 ~ (1 | fyear), data=sdat.lgr, family=binomial)
  tmp.sum.css2 = summary(tmp.fit.css2)
  # collect outputs
  simout.lgr.css2[k,] = c(as.vector(tmp.sum.css2$coef[2, -3]), AIC(tmp.fit.css2.yr), AIC(tmp.fit.css2))

  ## Bypass effect = Alternative value 1 
  # generate random number of returning fish
  set.seed(k*10)
  tmp.return.alt1 = rbinom(nobs.lgr, size=sdat.lgr$nfish, prob=sdat.lgr$psar.alt1)
  tmp.noreturn.alt1 = as.integer(sdat.lgr$nfish - tmp.return.alt1)
  tmp.adult_return.alt1 = cbind(tmp.return.alt1, tmp.noreturn.alt1)
  # fit models with and without nbypass
  tmp.fit.alt1 = glmer(tmp.adult_return.alt1 ~ nbypass + (1 | fyear), data=sdat.lgr, family=binomial)
  tmp.fit.alt1.yr = glmer(tmp.adult_return.alt1 ~ (1 | fyear), data=sdat.lgr, family=binomial)
  tmp.sum.alt1 = summary(tmp.fit.alt1)
  # collect outputs
  simout.lgr.alt1[k,] = c(as.vector(tmp.sum.alt1$coef[2, -3]), AIC(tmp.fit.alt1.yr), AIC(tmp.fit.alt1))

  ## Bypass effect = Alternative value 2 
  # generate random number of returning fish
  set.seed(k*10)
  tmp.return.alt2 = rbinom(nobs.lgr, size=sdat.lgr$nfish, prob=sdat.lgr$psar.alt2)
  tmp.noreturn.alt2 = as.integer(sdat.lgr$nfish - tmp.return.alt2)
  tmp.adult_return.alt2 = cbind(tmp.return.alt2, tmp.noreturn.alt2)
  # fit models with and without nbypass
  tmp.fit.alt2 = glmer(tmp.adult_return.alt2 ~ nbypass + (1 | fyear), data=sdat.lgr, family=binomial)
  tmp.fit.alt2.yr = glmer(tmp.adult_return.alt2 ~  (1 | fyear), data=sdat.lgr, family=binomial)
  tmp.sum.alt2 = summary(tmp.fit.alt2)
  # collect outputs
  simout.lgr.alt2[k,] = c(as.vector(tmp.sum.alt2$coef[2, -3]), AIC(tmp.fit.alt2.yr), AIC(tmp.fit.alt2))

}


## --- Calculate Power for Various Criteria ----
## -- alpha = 0.05
power.lgr.css1 = sum(simout.lgr.css1[,3]<=0.05)/nsim
power.lgr.css2 = sum(simout.lgr.css2[,3]<=0.05)/nsim
power.lgr.alt1 = sum(simout.lgr.alt1[,3]<=0.05)/nsim
power.lgr.alt2 = sum(simout.lgr.alt2[,3]<=0.05)/nsim

## -- alpha = 0.10
power.lgr.css1.10 = sum(simout.lgr.css1[,3]<=0.1)/nsim
power.lgr.css2.10 = sum(simout.lgr.css2[,3]<=0.1)/nsim
power.lgr.alt1.10 = sum(simout.lgr.alt1[,3]<=0.1)/nsim
power.lgr.alt2.10 = sum(simout.lgr.alt2[,3]<=0.1)/nsim

## -- delta AIC >= 2
power.lgr.css1.A2 = sum((simout.lgr.css1[,4]-simout.lgr.css1[,5])>=2)/nsim
power.lgr.css2.A2 = sum((simout.lgr.css2[,4]-simout.lgr.css2[,5])>=2)/nsim
power.lgr.alt1.A2 = sum((simout.lgr.alt1[,4]-simout.lgr.alt1[,5])>=2)/nsim
power.lgr.alt2.A2 = sum((simout.lgr.alt2[,4]-simout.lgr.alt2[,5])>=2)/nsim

## -- delta AIC >= 0
power.lgr.css1.A0 = sum((simout.lgr.css1[,4]-simout.lgr.css1[,5])>=0)/nsim
power.lgr.css2.A0 = sum((simout.lgr.css2[,4]-simout.lgr.css2[,5])>=0)/nsim
power.lgr.alt1.A0 = sum((simout.lgr.alt1[,4]-simout.lgr.alt1[,5])>=0)/nsim
power.lgr.alt2.A0 = sum((simout.lgr.alt2[,4]-simout.lgr.alt2[,5])>=0)/nsim



## Table for results
simsum.lgr = data.frame(species="Steelhead", tagsite="LGR", matrix(0,4,5))
names(simsum.lgr)[3:7] = c("parameter", "alpha_05", "alpha_10", "AIC2", "AIC0")
simsum.lgr$parameter = c(tmp.nbyp.css1,tmp.nbyp.css2, tmp.nbyp.alt1,tmp.nbyp.alt2)
simsum.lgr$alpha_05 = c(power.lgr.css1, power.lgr.css2, power.lgr.alt1, power.lgr.alt2)
simsum.lgr$alpha_10 = c(power.lgr.css1.10, power.lgr.css2.10, power.lgr.alt1.10, power.lgr.alt2.10)
simsum.lgr$AIC2 = c(power.lgr.css1.A2, power.lgr.css2.A2, power.lgr.alt1.A2, power.lgr.alt2.A2)
simsum.lgr$AIC0 = c(power.lgr.css1.A0, power.lgr.css2.A0, power.lgr.alt1.A0, power.lgr.alt2.A0)
simsum.lgr = simsum.lgr[order(simsum.lgr$parameter, decreasing=T), ]


##  Output Combined Results
simsum.out = rbind(simsum.ulgr, simsum.lgr)
write.csv(simsum.out, file="output/steelhead_powersim_summary_table.csv", row.names=F, quote=F)

 


## Try with fixed effects for year
# NOTE: results for this were referenced in the response to comment but not directly reported
# This was tested as an explanation for lower power reported by Storch et al
nsim = 1000
simout.lgr.css1.fixef = matrix(0, nsim, 3)

for (k in 1:nsim) {
  cat("running sim", k, "of", nsim, "\n")
  
  set.seed(k*10)
  tmp.return.css = rbinom(nobs.lgr, size=sdat.lgr$nfish, prob=sdat.lgr$psar.css1)
  tmp.noreturn.css = as.integer(sdat.lgr$nfish - tmp.return.css)
  tmp.adult_return.css = cbind(tmp.return.css, tmp.noreturn.css)

  tmp.fit.css = glm(tmp.adult_return.css ~ nbypass + fyear, data=sdat.lgr, family=binomial)
  tmp.sum.css = summary(tmp.fit.css)
  simout.lgr.css1.fixef[k,] = as.vector(tmp.sum.css$coef[2, -3])

}

power.lgr.css1.fixef = sum(simout.lgr.css1.fixef[,3]<=0.05)/nsim
power.lgr.css1.fixef




### -----------------------------------------------------
##      Power by Sample Size
### ------------------------------------------------------

# The following simulation results were used for power curves in Figure 2


## ----- Get Expected cell probs ------

# function for multinomial negative log likelihood
# assumes same detection across all dams in a year
multiNLL <- function(parm, ncell, counts) {
  pval <- plogis(parm)
  pvec <- dbinom(0:ncell, size=ncell, prob=pval)
  llik <- dmultinom(counts, size=sum(counts), prob=pvec, log=T)
  return(-llik)
}

## Estimate shared detection probabilities by year
pinit = qlogis(0.5)
pout.ulgr = numeric(nyear)
for (i in 1:nyear) {
  tmp.pfit = optim(par=pinit, fn=multiNLL, method="BFGS", ncell=7, counts=byptab.ulgr[i,])
  pout.ulgr[i] = plogis(tmp.pfit$par)
}

## Use detection probs to generate probabilities for 0 to 7 bypass events
pmat.ulgr = matrix(0, nyear, 8)
for (i in 1:nyear) {
  pmat.ulgr[i, ] = dbinom(0:7, size=7, prob=pout.ulgr[i])
}

## baseline number of fish in each year
nvec.ulgr = as.vector(apply(byptab.ulgr, 1, sum))


### ---- Run Simulations -----------
nmult = 1:9  # multiplier for sample sizes
nsim = 500   # number of simulations at each sample size level
npowvec.ulgr = numeric(length(nmult))  #vector to hold results

for (j in nmult) {
 
  tmp.simout = matrix(0, nsim, 3)  #parm, se, pval, aic.yronly, aic.full

  for (k in 1:nsim) {
    cat("mult = ",j,"sim =", k, "of", nsim, "\n")

    # Generate number of bypass events in each year
    set.seed(k*9)
    tmp.sdat1.ulgr = data.frame(year=as.character(rep(yrset, rep(8,nyear))), nbypass=rep(0:7,nyear), nfish=0)
    tmp.nfv.ulgr = numeric()
    for (i in 1:nyear) {
      tmp.nf = rmultinom(1, size=nmult[j]*nvec.ulgr[i], prob=pmat.ulgr[i,])
      tmp.nfv.ulgr = c(tmp.nfv.ulgr, tmp.nf)
    }
    # Generate return probabilities
    tmp.sdat1.ulgr$nfish = tmp.nfv.ulgr
    tmp.sdat1.ulgr$nu.css1 = tmp.int.ulgr + tmp.yreff.ulgr + tmp.nbyp.css1*tmp.sdat1.ulgr$nbypass
    tmp.sdat1.ulgr$psar.css1 = plogis(tmp.sdat1.ulgr$nu.css1)
    tmp.sdat.ulgr = tmp.sdat1.ulgr[tmp.sdat1.ulgr$nfish!=0, ]
    tmp.sdat.ulgr$fyear = as.character(tmp.sdat.ulgr$year)
    tmp.nobs.ulgr = nrow(tmp.sdat.ulgr)
    # Generate number of fish returning  
    set.seed(k*10)
    tmp.return.css1 = rbinom(tmp.nobs.ulgr, size=tmp.sdat.ulgr$nfish, prob=tmp.sdat.ulgr$psar.css1)
    tmp.noreturn.css1 = as.integer(tmp.sdat.ulgr$nfish - tmp.return.css1)
    tmp.adult_return.css1 = cbind(tmp.return.css1, tmp.noreturn.css1)
    # Fit model and record pvalue for nbypass
    tmp.fit.css1 = glmer(tmp.adult_return.css1 ~ nbypass + (1 | fyear), data=tmp.sdat.ulgr, family=binomial)
    tmp.sum.css1 = summary(tmp.fit.css1)
    tmp.simout[k,] = as.vector(tmp.sum.css1$coef[2, -3]) 
  }

 # Calculate power
 npowvec.ulgr[j] = sum(tmp.simout[,3]<=0.05)/nsim

}


## Write outputs to file for plot generation
write.csv(data.frame(size=1:9, power=npowvec.ulgr), "output/steelhead_power_by_samplesize_summary.csv", quote=F, row.names=F)







