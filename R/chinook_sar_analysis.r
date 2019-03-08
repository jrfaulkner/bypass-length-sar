## Chinook SAR Analysis
## For model fitting and simulations in support of:
## Faulkner, JR, BL Bellerud, DL Widener, and RW Zabel. 2019. Associations among fish length, dam passage history, and survival to adulthood in two at-risk species of Pacific salmon.



# These packages need to be installed and loaded
library(lme4)
library(Matrix)
library(mvtnorm)




###-------------------------------------------------------------------------------------
##                 Data
###-------------------------------------------------------------------------------------

### --- Read data files
cdat.ulgr <- read.csv("../../data/chinook_ULGR_SAR_data.csv", stringsAsFactors=F)
cdat.lgr <- read.csv("../../data/chinook_LGR_SAR_data.csv", stringsAsFactors=F)


## --- Make new variables
## make year a factor
cdat.ulgr$fyear <- as.factor(cdat.ulgr$year)
cdat.lgr$fyear <- as.factor(cdat.lgr$year)



###----------------------------------------------------------------------------------------
##            Model Fitting
###----------------------------------------------------------------------------------------


### ----  Functions ---------------------------------

## Function to create table of AIC ranks and other stats
getAICsummary <- function(fit.list, mod.list, modnums=NULL) {
  nmod <- length(fit.list)-1
  if (is.null(modnums)) modnums <- 1:nmod
  aicv1 <- as.vector((sapply(fit.list[1:nmod], FUN=extractAIC))[2,])  #extract aic
  minaic1 <- min(aicv1)	#min aic
  delta1 <- aicv1 - minaic1	#delta aic
  mod.lik1 <- exp(-0.5*delta1)	#likelihood of model
  aicw1 <- mod.lik1/sum(mod.lik1)	#aic weight
  evratio1 <- 1/mod.lik1		# evidence ratio
  aicm1a <- cbind(modnums, aicv1, aicw1)
  aicm1a <- aicm1a[order(aicm1a[,2]), ]  #sort by aic
  aicm1b <- cbind(aicm1a, cumsum(aicm1a[,3]), 1:nmod)	#assign rank
  aicm1b <- aicm1b[order(aicm1b[,1]), ]      #resort by model #
  aicout1 <- data.frame(aicm1b[,c(1:2,5)] , delta1, aicm1b[,3:4], mod.lik1, evratio1, mod.list, stringsAsFactors=F)  #combine all into data frame
  names(aicout1) <- c("Model", "AIC", "Rank", "DeltaAIC", "Weight","CumWeight", "ModelLik", "EvRatio","ModelForm")
  aicout1 <- aicout1[order(aicout1$Rank), ]  #resort by rank
  return(aicout1)
}


# function for extracting odds ratio and CI for single fixed effect from logistic regression
extractCI <- function(modsum, varname, scale="odds", alpha=0.05){
	vals <- modsum$coef[varname, 1:2]
	est <- vals[1]
	se <- vals[2]
	ci <- est +c(-1,1)*qnorm(p=1-alpha/2)*se
	if (scale=="odds") out <- c(estimate=exp(est), ci.low=exp(ci[1]), ci.up=exp(ci[2]) )
	if (scale=="log.odds") out <- c(estimate=est, ci.low=ci[1], ci.up=ci[2] )
	return(out)
}




#### -------  Model Formulas ---------------------------

### Variable shorthand
# Y = random effect for migration year
# C = reartype (clipped or not)
# D = last day of detection
# D2 = D squared
# R = release site (for ULGR only)
# S = site of last detection (BON or TWX)
# L = length
# B = bypass yes/no
# N = number of bypasses
# NC = categorical number of bypasses

### Possible covariate-only models
# note: all models have Y and C

# ULGR models
# 1) Y + C
# 2) Y + C + D
# 3) Y + C + D + D2
# 4) Y + C + R
# 5) Y + C + S
# 6) Y + C + D + R
# 7) Y + C + D + S
# 8) Y + C + D + D2 + R
# 9) Y + C + D + D2 + S
# 10) Y + C + R + S
# 11) Y + C + D + R + S
# 12) Y + C + D + D2 + R + S
# 13) Y + C + D + Y*D
# 14) Y + C + D + D2 + Y*D
# 15) Y + C + D + D2 + Y*D + Y*D2
# 16) Y + C + D + R + Y*D
# 17) Y + C + D + S + Y*D
# 18) Y + C + D + D2 + R + Y*D 
# 19) Y + C + D + D2 + R + Y*D + Y*D2
# 20) Y + C + D + D2 + S + Y*D 
# 21) Y + C + D + D2 + S + Y*D + Y*D2
# 22) Y + C + D + R + S + Y*D
# 23) Y + C + D + D2 + R + S + Y*D 
# 24) Y + C + D + D2 + R + S + Y*D + Y*D2


# LGR models
# 1) Y + C
# 2) Y + C + D
# 3) Y + C + D + D2
# 4) Y + C + S
# 5) Y + C + D + S
# 6) Y + C + D + D2 + S
# 7) Y + C + D + Y*D
# 8) Y + C + D + D2 + Y*D
# 9) Y + C + D + D2 + Y*D + Y*D2
# 10) Y + C + D + S + Y*D
# 11) Y + C + D + D2 + S + Y*D 
# 12) Y + C + D + D2 + S + Y*D + Y*D2


## models with length and/or bypass variables
# cov + L
# cov + B
# cov + N
# cov + NC
# cov + L + B
# cov + L + N
# cov + L + NC


covsets.ulgr <- c("clipped + (1 | fyear)",
  "clipped + zlast_day + (1 | fyear)",   
  "clipped + zlast_day + zlast_day_sq + (1 | fyear)",
  "clipped + rel_site2 + (1 | fyear)",
  "clipped + last_site + (1 | fyear)", 
  "clipped + zlast_day + rel_site2 + (1 | fyear)",
  "clipped + zlast_day + last_site + (1 | fyear)",
  "clipped + zlast_day + zlast_day_sq + rel_site2 + (1 | fyear)",
  "clipped + zlast_day + zlast_day_sq + last_site + (1 | fyear)",
  "clipped + rel_site2 + last_site + (1 | fyear)", 
  "clipped + zlast_day + rel_site2 + last_site + (1 | fyear)", 
  "clipped + zlast_day + zlast_day_sq + rel_site2 + last_site + (1 | fyear)", 
  "clipped + zlast_day + (zlast_day | fyear)",   
  "clipped + zlast_day + zlast_day_sq + (zlast_day | fyear)",   
  "clipped + zlast_day + zlast_day_sq + (zlast_day + zlast_day_sq | fyear)",   
  "clipped + zlast_day + rel_site2 + (zlast_day | fyear)",   
  "clipped + zlast_day + last_site + (zlast_day | fyear)",   
  "clipped + zlast_day + zlast_day_sq + rel_site2 + (zlast_day | fyear)",   
  "clipped + zlast_day + zlast_day_sq + rel_site2 + (zlast_day + zlast_day_sq | fyear)",   
  "clipped + zlast_day + zlast_day_sq + last_site + (zlast_day | fyear)",   
  "clipped + zlast_day + zlast_day_sq + last_site + (zlast_day + zlast_day_sq | fyear)",   
  "clipped + zlast_day + rel_site2 + last_site + (zlast_day | fyear)", 
  "clipped + zlast_day + zlast_day_sq + rel_site2 + last_site + (zlast_day | fyear)", 
  "clipped + zlast_day + zlast_day_sq + rel_site2 + last_site + (zlast_day + zlast_day_sq  | fyear)")


covsets.lgr <- c("clipped + (1 | fyear)",
  "clipped + zlast_day + (1 | fyear)",   
  "clipped + zlast_day + zlast_day_sq + (1 | fyear)",
  "clipped + last_site + (1 | fyear)", 
  "clipped + zlast_day + last_site + (1 | fyear)",
  "clipped + zlast_day + zlast_day_sq + last_site + (1 | fyear)",
  "clipped + zlast_day + (zlast_day | fyear)",   
  "clipped + zlast_day + zlast_day_sq + (zlast_day | fyear)",   
  "clipped + zlast_day + zlast_day_sq + (zlast_day + zlast_day_sq | fyear)",   
  "clipped + zlast_day + last_site + (zlast_day | fyear)",   
  "clipped + zlast_day + zlast_day_sq + last_site + (zlast_day | fyear)",   
  "clipped + zlast_day + zlast_day_sq + last_site + (zlast_day + zlast_day_sq | fyear)") 


testsets <- c("zlength", "bypassed", "nbypass_comb4", "nbypass", "zlength + bypassed", "zlength + nbypass_comb4", "zlength + nbypass")

covformulas.ulgr <- covformulas.lgr <- list()
nmods.ulgr <- length(covsets.ulgr)
nmods.lgr <- length(covsets.lgr)

for (i in 1:nmods.ulgr) {
	covformulas.ulgr[[i]] <- formula(paste("adult_return ~ ", covsets.ulgr[i] ))
}

for (i in 1:nmods.lgr) {
	covformulas.lgr[[i]] <- formula(paste("adult_return ~ ", covsets.lgr[i] ))
}




###-------------------------------------------------
##     Released Upstream of Lower Granite
###-------------------------------------------------


### -------   Fit Covariate-Only Models ------------

modfits.cov.ulgr <- list()
convcheck.cov.ulgr <- numeric(nmods.ulgr)

for (k in 1:nmods.ulgr) {
	cat("fitting model: ", k, " of ", nmods.ulgr, "\n" )
	tmp.fit <- glmer(covformulas.ulgr[[k]], family=binomial, data=cdat.ulgr)
  if (length(tmp.fit@optinfo$conv$lme4)==0) convcheck.cov.ulgr[k] <- 0
  if (length(tmp.fit@optinfo$conv$lme4)!=0) convcheck.cov.ulgr[k] <- tmp.fit@optinfo$conv$lme4$code
  modfits.cov.ulgr[[k]] <- tmp.fit 
	# if did not converge, try to refit
	if (convcheck.cov.ulgr[k]!=0) {
		cat("re-fitting model: ", k, " of ", nmods.ulgr, "\n" )
		tmp.ss <- getME(tmp.fit,c("theta","fixef"))
		tmp.fit2 <- update(tmp.fit,start=tmp.ss,control=glmerControl(optCtrl=list(maxfun=2e4)))
	  if (length(tmp.fit2@optinfo$conv$lme4)==0) convcheck.cov.ulgr[k] <- 0
	  if (length(tmp.fit2@optinfo$conv$lme4)!=0) convcheck.cov.ulgr[k] <- tmp.fit2@optinfo$conv$lme4$code
	  modfits.cov.ulgr[[k]] <- tmp.fit2 
	}
}


# check if any models did not converge (any non-zero)
convcheck.cov.ulgr
# all are fine

# Calculate AIC
aicsum.cov.ulgr <- getAICsummary(modfits.cov.ulgr, covsets.ulgr)

# find best model
bestmod.cov.ulgr <- aicsum.cov.ulgr$Model[1]


## --------- Fit length and bypass models ------------------------

# Make model formulas for bypass and length models
formulas.ulgr <- list()
nmods <- length(testsets)
for (i in 1:nmods) {
  formulas.ulgr[[i]] <- 	formula(paste("adult_return ~ ", testsets[i],  "+", covsets.ulgr[bestmod.cov.ulgr] ))
}

modfits.ulgr <- list()
modfits.ulgr[[1]] <- modfits.cov.ulgr[[bestmod.cov.ulgr]]
convcheck.ulgr <- numeric(nmods+1)

for (k in 1:nmods) {
	cat("fitting model: ", k, " of ", nmods, "\n" )
	tmp.fit <- glmer(formulas.ulgr[[k]], family=binomial, data=cdat.ulgr)
  if (length(tmp.fit@optinfo$conv$lme4)==0) convcheck.ulgr[k+1] <- 0
  if (length(tmp.fit@optinfo$conv$lme4)!=0) convcheck.ulgr[k+1] <- tmp.fit@optinfo$conv$lme4$code
  modfits.ulgr[[k+1]] <- tmp.fit 
	# if did not converge, try to refit
	if (convcheck.ulgr[k+1]!=0) {
		cat("re-fitting model: ", k, " of ", nmods, "\n" )
		tmp.ss <- getME(tmp.fit,c("theta","fixef"))
		tmp.fit2 <- update(tmp.fit,start=tmp.ss,control=glmerControl(optCtrl=list(maxfun=2e4)))
	  if (length(tmp.fit2@optinfo$conv$lme4)==0) convcheck.ulgr[k+1] <- 0
	  if (length(tmp.fit2@optinfo$conv$lme4)!=0) convcheck.ulgr[k+1] <- tmp.fit2@optinfo$conv$lme4$code
	  modfits.ulgr[[k+1]] <- tmp.fit2 
	}
}


# check if any models did not converge (any non-zero)
convcheck.ulgr
# all are fine

# Calculate AIC
aicsum.ulgr <- getAICsummary(modfits.ulgr, c("cov", paste("cov +", testsets ) ) )
aicsum.ulgr


# get model summaries
modsums.ulgr <- lapply(modfits.ulgr, summary)

# make a table of relvant pvalues
pvals.ulgr <- data.frame(Model=1:(nmods+1), matrix(NA, nrow=(nmods+1), ncol=2) )
names(pvals.ulgr)[2:3] <- c("p.length", "p.bypass")
for (k in 2:(nmods+1) ) {
	if (k %in% c(2,6:8)) pvals.ulgr[k, 2] <- round(modsums.ulgr[[k]]$coef["zlength", 4], 3) 
	if (k %in% c(3,6))   pvals.ulgr[k, 3] <- round(modsums.ulgr[[k]]$coef["bypassed", 4], 3) 
	if (k %in% c(5,8))   pvals.ulgr[k, 3] <- round(modsums.ulgr[[k]]$coef["nbypass", 4], 3) 
	if (k==4) pvals.ulgr[k, 3] <- round(anova(modfits.ulgr[[4]], modfits.ulgr[[1]])$P[2], 3)
	if (k==7) pvals.ulgr[k, 3] <- round(anova(modfits.ulgr[[7]], modfits.ulgr[[2]])$P[2], 3)
}


# Create summary stats for Table 2
tab2.ulgr <- merge(aicsum.ulgr[order(aicsum.ulgr$Model), c("Model", "ModelForm", "DeltaAIC", "Weight") ], pvals.ulgr)
tab2.ulgr[, 3] <- round(tab2.ulgr[, 3], 1)
tab2.ulgr[, 4] <- round(tab2.ulgr[, 4], 2)

write.csv(tab2.ulgr, "output/chinook_ulgr_table2.csv", row.names=F, quote=F)


# --- Get CI's for odds ratios on length and nbypass

# zlength
extractCI(modsums.ulgr[[8]], "zlength", scale="log.odds")

# nbypass
extractCI(modsums.ulgr[[8]], "nbypass", scale="log.odds")





###-------------------------------------------------
##     Released at Lower Granite
###-------------------------------------------------


### -------   Fit Covariate-Only Models ------------

modfits.cov.lgr <- list()
convcheck.cov.lgr <- numeric(nmods.lgr)

for (k in 1:nmods.lgr) {
	cat("fitting model: ", k, " of ", nmods.lgr, "\n" )
	tmp.fit <- glmer(covformulas.lgr[[k]], family=binomial, data=cdat.lgr)
  if (length(tmp.fit@optinfo$conv$lme4)==0) convcheck.cov.lgr[k] <- 0
  if (length(tmp.fit@optinfo$conv$lme4)!=0) convcheck.cov.lgr[k] <- tmp.fit@optinfo$conv$lme4$code
  modfits.cov.lgr[[k]] <- tmp.fit 
	# if did not converge, try to refit
	if (convcheck.cov.lgr[k]!=0) {
		cat("re-fitting model: ", k, " of ", nmods.lgr, "\n" )
		tmp.ss <- getME(tmp.fit,c("theta","fixef"))
		tmp.fit2 <- update(tmp.fit,start=tmp.ss,control=glmerControl(optCtrl=list(maxfun=2e4)))
	  if (length(tmp.fit2@optinfo$conv$lme4)==0) convcheck.cov.lgr[k] <- 0
	  if (length(tmp.fit2@optinfo$conv$lme4)!=0) convcheck.cov.lgr[k] <- tmp.fit2@optinfo$conv$lme4$code
	  modfits.cov.lgr[[k]] <- tmp.fit2 
	}
}


# check if any models did not converge (any non-zero)
convcheck.cov.lgr
# all are fine

# Calculate AIC
aicsum.cov.lgr <- getAICsummary(modfits.cov.lgr, covsets.lgr)

# find best model
bestmod.cov.lgr <- aicsum.cov.lgr$Model[1]



## --------- Fit length and bypass models ------------------------

# Make model formulas for bypass and length models
formulas.lgr <- list()
nmods <- length(testsets)
for (i in 1:nmods) {
  formulas.lgr[[i]] <- 	formula(paste("adult_return ~ ", testsets[i],  "+", covsets.lgr[bestmod.cov.lgr] ))
}

modfits.lgr <- list()
modfits.lgr[[1]] <- modfits.cov.lgr[[bestmod.cov.lgr]]
convcheck.lgr <- numeric(nmods+1)

for (k in 1:nmods) {
	cat("fitting model: ", k, " of ", nmods, "\n" )
	tmp.fit <- glmer(formulas.lgr[[k]], family=binomial, data=cdat.lgr)
  if (length(tmp.fit@optinfo$conv$lme4)==0) convcheck.lgr[k+1] <- 0
  if (length(tmp.fit@optinfo$conv$lme4)!=0) convcheck.lgr[k+1] <- tmp.fit@optinfo$conv$lme4$code
  modfits.lgr[[k+1]] <- tmp.fit 
	# if did not converge, try to refit
	if (convcheck.lgr[k+1]!=0) {
		cat("re-fitting model: ", k, " of ", nmods, "\n" )
		tmp.ss <- getME(tmp.fit,c("theta","fixef"))
		tmp.fit2 <- update(tmp.fit,start=tmp.ss,control=glmerControl(optCtrl=list(maxfun=2e4)))
	  if (length(tmp.fit2@optinfo$conv$lme4)==0) convcheck.lgr[k+1] <- 0
	  if (length(tmp.fit2@optinfo$conv$lme4)!=0) convcheck.lgr[k+1] <- tmp.fit2@optinfo$conv$lme4$code
	  modfits.lgr[[k+1]] <- tmp.fit2 
	}
}


# check if any models did not converge (any non-zero)
convcheck.lgr
# all are fine

# Calculate AIC
aicsum.lgr <- getAICsummary(modfits.lgr, c("cov", paste("cov +", testsets ) ) )
aicsum.lgr


# get model summaries
modsums.lgr <- lapply(modfits.lgr, summary)

# make a table of relvant pvalues
pvals.lgr <- data.frame(Model=1:(nmods+1), matrix(NA, nrow=(nmods+1), ncol=2) )
names(pvals.lgr)[2:3] <- c("p.length", "p.bypass")
for (k in 2:(nmods+1) ) {
	if (k %in% c(2,6:8)) pvals.lgr[k, 2] <- round(modsums.lgr[[k]]$coef["zlength", 4], 3) 
	if (k %in% c(3,6))   pvals.lgr[k, 3] <- round(modsums.lgr[[k]]$coef["bypassed", 4], 3) 
	if (k %in% c(5,8))   pvals.lgr[k, 3] <- round(modsums.lgr[[k]]$coef["nbypass", 4], 3) 
	if (k==4) pvals.lgr[k, 3] <- round(anova(modfits.lgr[[4]], modfits.lgr[[1]])$P[2], 3)
	if (k==7) pvals.lgr[k, 3] <- round(anova(modfits.lgr[[7]], modfits.lgr[[2]])$P[2], 3)
}


# Create summary stats for Table 2
tab2.lgr <- merge(aicsum.lgr[order(aicsum.lgr$Model), c("Model", "ModelForm", "DeltaAIC", "Weight") ], pvals.lgr)
tab2.lgr[, 3] <- round(tab2.lgr[, 3], 1)
tab2.lgr[, 4] <- round(tab2.lgr[, 4], 2)

write.csv(tab2.lgr, "output/chinook_lgr_table2.csv", row.names=F, quote=F)


# --- Get CI's for odds ratios on length and nbypass

# zlength
extractCI(modsums.lgr[[8]], "zlength", scale="log.odds")

# nbypass
extractCI(modsums.lgr[[8]], "nbypass", scale="log.odds")





####---------------------------------------------------------------------------------------------------
##                       Simulations
####---------------------------------------------------------------------------------------------------


## Objective is to generate data from specific models and fit different models to that data to test bypass effect.
## Models of interest (as labeled in the manuscript):
# M0: covs only
# M1: covs + nbypass
# M2: covs + nbypass + length
# M3: covs + length

## Three comparisons of interest are:
# 1) Generate from M0 and fit M1, test nbypass effect
# 2) Generate from M3 and fit M1, test nbypass effect
# 3) Generate from M3 and fit M2, test nbypass effect


# NOTE: M0 = best covariate model
#       M1 = model 4 in model formulas and model 5 in modelfits
#       M2 = model 7 in model formulas and model 8 in modelfits
#       M3 = model 1 in model formulas and model 2 in modelfits


## !!! NOTE:  These simulations can take many days to run completely.


### ----  Functions ---------------------------------

fullCov <- function(fitobj, gmat) {
  fpred <- predict(fitobj, type="response")
  fX <- as(getME(fitobj, "X"), Class="Matrix")						
  fZ <- as(getME(fitobj, "Z"), Class="Matrix")		
  fZt <- t(fZ)		
  fXt <- t(fX)
  fSi <- Diagonal(x=fpred*(1-fpred) )
  fXtSi <- crossprod(fX, fSi)
  fZtSi <- crossprod(fZ, fSi)
  fG <- gmat
  cvout <- as(solve(rBind( cBind(tcrossprod(fXtSi, fXt), tcrossprod(fXtSi, fZt)  ),
cBind(tcrossprod(fZtSi, fXt), tcrossprod(fZtSi, fZt) + solve(fG) ) )  ), Class="matrix")
  return(cvout)
}

interleave <- function(v1,v2){
  ord1 <- 2*(1:length(v1))-1
  ord2 <- 2*(1:length(v2))
  c(v1,v2)[order(c(ord1,ord2))]
}

## Function to extract conditional rand cov matrix
condVar <- function(object) {
  s2 <- sigma(object)^2
  Lamt <- getME(object, "Lambdat")
  L <- getME(object, "L")
  LL <- solve(L, Lamt, system = "A")
  s2 * crossprod(Lamt, LL)
}







####---------------------------------------------------------
##     Tagged Upstream of Lower Granite
####---------------------------------------------------------


###-----------------------------
##     Generate Data from M0
###-----------------------------


## ---- Set up covariance matrices and parameter means------
G.ulgr.m0 <- Diagonal(n=11, x=VarCorr(modfits.ulgr[[1]])$fyear[1]) 
X.ulgr.m0 <- as(getME(modfits.ulgr[[1]], "X"), Class="Matrix")						
Z.ulgr.m0 <- as(getME(modfits.ulgr[[1]], "Z"), Class="Matrix")						 
VC.ulgr.m0 <- fullCov(modfits.ulgr[[1]], G.ulgr.m0)
bmu.ulgr.m0 <- fixef(modfits.ulgr[[1]])
rand.ulgr.m0 <- as.matrix(ranef(modfits.ulgr[[1]])$fyear)
rmu.ulgr.m0 <- as.vector(rand.ulgr.m0)


### --- Run simulations ------------------------

formula.ulgr.M0 <- formula(paste("adult_return ~ ", covsets.ulgr[[bestmod.cov.ulgr]] ))

nsim <- 1200

simout.ulgr.m0 <- data.frame(matrix(NA, nsim, 7 ))
names(simout.ulgr.m0) <- c("est.nbyp.M1", "SE.nbyp.M1", "zval.nbyp.M1", "pval.nbyp.M1","aic.M1", "aic.M0", "converge")
simout.ulgr.m0$converge=TRUE

tmpdat.ulgr <- cdat.ulgr

for (k in 1:nsim) {
  set.seed(k*10) 
  cat("sim ", k, " of", nsim, "\n")
  tmp.pars <- rmvnorm(n=1, mean=c(bmu.ulgr.m0, rmu.ulgr.m0), sigma=VC.ulgr.m0 )
  tmp.pred.eta <- X.ulgr.m0%*%as(as.matrix(tmp.pars[1:5]), Class="Matrix") + Z.ulgr.m0%*%as(as.matrix(tmp.pars[6:16]), Class="Matrix")
  tmp.pred <- plogis(as.vector(as(tmp.pred.eta, Class="matrix")))
  # generate new adult return data
  tmp.y <- rbinom(n=length(tmp.pred.eta), size=1, prob=tmp.pred)
  tmpdat.ulgr$adult_return <- tmp.y

  # fit models
  cat("\tfitting model M0", "\n")
  tmpfit.M0 <- glmer(formula.ulgr.M0, family=binomial, data=tmpdat.ulgr)
  cat("\tfitting model M1", "\n")
  tmpfit.M1 <- glmer(formulas.ulgr[[4]], family=binomial, data=tmpdat.ulgr)

  # re-fit any that did not converge
  if (length(tmpfit.M0@optinfo$conv$lme4)!=0){
    cat("\tre-fitting model M0", "\n")
    tmp.ss.M0 <- getME(tmpfit.M0,c("theta","fixef"))
    tmpfit.M0 <- update(tmpfit.M0, start=tmp.ss.M0,control=glmerControl(optCtrl=list(maxfun=2e4)))
  }
  if (length(tmpfit.M1@optinfo$conv$lme4)!=0){
    cat("\tre-fitting model M1", "\n")
    tmp.ss.M1 <- getME(tmpfit.M1,c("theta","fixef"))
    tmpfit.M1 <- update(tmpfit.M1, start=tmp.ss.M1,control=glmerControl(optCtrl=list(maxfun=2e4)))
  }
	
  # if any still not converged then mark this sim as not converged
  if (any(!is.null(tmpfit.M0@optinfo$conv$lme4$messages), !is.null(tmpfit.M1@optinfo$conv$lme4$messages)) ) {
    simout.ulgr.m0$converge[k]=FALSE		
  }

  # collect outputs
  simout.ulgr.m0[k, 1:4] <- summary(tmpfit.M1)$coef["nbypass", ]
  simout.ulgr.m0[k, 5] <- AIC(tmpfit.M1)
  simout.ulgr.m0[k, 6] <- AIC(tmpfit.M0)

  write.csv(simout.ulgr.m0, file="output/chinook_ulgr_dataM0_simulation_outputs.csv", quote=FALSE, row.names=FALSE)
  if (k==nsim) cat("Done!", "\n")
}


### ----- Summarize outputs -------------------------

## remove any non-converged runs 
simout.ulgr.m0f <- simout.ulgr.m0
simout.ulgr.m0 <- simout.ulgr.m0[simout.ulgr.m0$converge==TRUE & !is.na(simout.ulgr.m0$est.nbyp.M1), ]
# check if there are at least 1000 remaining sims first
dim(simout.ulgr.m0)

## reduce to first 1000 of remaining and run additional if necessary
nsim <- 1000
simout.ulgr.m0 <- simout.ulgr.m0[1:nsim, ]

## create outputs for Table 3
p.n.ulg.m0.m1 <- round(100*sum(simout.ulgr.m0$est.nbyp.M1 < 0 )/nsim, 1)
p.ns.ulg.m0.m1 <- round(100*sum(simout.ulgr.m0$est.nbyp.M1 < 0 & simout.ulgr.m0$pval.nbyp.M1 <= 0.05 )/nsim, 1)
p.na.ulg.m0.m1 <- round(100*sum(simout.ulgr.m0$est.nbyp.M1 < 0 & simout.ulgr.m0$aic.M1 < simout.ulgr.m0$aic.M0)/nsim, 1)

outcome.labs <- c("nbyp.neg", "nbyp.neg & p<0.05", "nbyp.neg & lower AIC")
expected.labs <- c(50.0, 2.5, 7.8)
tab3.ulgr.m0 <- data.frame(data.model=rep("M0", 3), fit.model=rep("M1",3), 
                    outcome=outcome.labs, expected=expected.labs,
                    chinook.ULGD=c(p.n.ulg.m0.m1, p.ns.ulg.m0.m1, p.na.ulg.m0.m1) )





###----------------------------------
##      Generate Data from M3
###----------------------------------


## ---- Set up covariance matrices and parameter means------
rcor.ulgr.m3 <-  as.data.frame(VarCorr(modfits.ulgr[[2]]))
G.ulgr.m3 <- Diagonal(n=11, x=VarCorr(modfits.ulgr[[2]])$fyear[1]) 
X.ulgr.m3 <- as(getME(modfits.ulgr[[2]], "X"), Class="Matrix")						
Z.ulgr.m3 <- as(getME(modfits.ulgr[[2]], "Z"), Class="Matrix")						 
VC.ulgr.m3 <- fullCov(modfits.ulgr[[2]], G.ulgr.m3)
bmu.ulgr.m3 <- fixef(modfits.ulgr[[2]])
rand.ulgr.m3 <- as.matrix(ranef(modfits.ulgr[[2]])$fyear)
rmu.ulgr.m3 <- as.vector(rand.ulgr.m3)



### --- Run simulations ------------------------

formula.ulgr.M0 <- formula(paste("adult_return ~ ", covsets.ulgr[[bestmod.cov.ulgr]] ))

nsim <- 1200  

simout.ulgr.m3 <- data.frame(matrix(NA, nsim, 13 ))
names(simout.ulgr.m3) <- c("est.nbyp.M1", "SE.nbyp.M1", "zval.nbyp.M1", "pval.nbyp.M1","aic.M1","est.nbyp.M2", "SE.nbyp.M2", "zval.nbyp.M2", "pval.nbyp.M2","aic.M2","aic.M3", "aic.M0",  "converge")
simout.ulgr.m3$converge=TRUE

tmpdat.ulgr <- cdat.ulgr

for (k in 1:nsim) {
	set.seed(k*10) 
  cat("sim ", k, " of", nsim, "\n")
  tmp.pars <- rmvnorm(n=1, mean=c(bmu.ulgr.m3, rmu.ulgr.m3), sigma=VC.ulgr.m3 )
  tmp.pred.eta <- X.ulgr.m3%*%as(as.matrix(tmp.pars[1:6]), Class="Matrix") + Z.ulgr.m3%*%as(as.matrix(tmp.pars[7:17]), Class="Matrix")
  tmp.pred <- plogis(as.vector(as(tmp.pred.eta, Class="matrix")))
  # generate new adult return data
  tmp.y <- rbinom(n=length(tmp.pred.eta), size=1, prob=tmp.pred)
  tmpdat.ulgr$adult_return <- tmp.y

  # fit models
  cat("\tfitting model M0", "\n")
  tmpfit.M0 <- glmer(formula.ulgr.M0, family=binomial, data=tmpdat.ulgr)
  cat("\tfitting model M1", "\n")
  tmpfit.M1 <- glmer(formulas.ulgr[[4]], family=binomial, data=tmpdat.ulgr)
  cat("\tfitting model M2", "\n")
  tmpfit.M2 <- glmer(formulas.ulgr[[7]], family=binomial, data=tmpdat.ulgr)
  cat("\tfitting model M3", "\n")
  tmpfit.M3 <- glmer(formulas.ulgr[[1]], family=binomial, data=tmpdat.ulgr)

  # re-fit any that did not converge
	if (length(tmpfit.M0@optinfo$conv$lme4)!=0){
	  cat("\tre-fitting model M0", "\n")
		tmp.ss.M0 <- getME(tmpfit.M0,c("theta","fixef"))
		tmpfit.M0 <- update(tmpfit.M0, start=tmp.ss.M0,control=glmerControl(optCtrl=list(maxfun=2e4)))
	}
	if (length(tmpfit.M1@optinfo$conv$lme4)!=0){
	  cat("\tre-fitting model M1", "\n")
		tmp.ss.M1 <- getME(tmpfit.M1,c("theta","fixef"))
		tmpfit.M1 <- update(tmpfit.M1, start=tmp.ss.M1,control=glmerControl(optCtrl=list(maxfun=2e4)))
	}
	if (length(tmpfit.M2@optinfo$conv$lme4)!=0){
	  cat("\tre-fitting model M2", "\n")
		tmp.ss.M2 <- getME(tmpfit.M2,c("theta","fixef"))
		tmpfit.M2 <- update(tmpfit.M2, start=tmp.ss.M2, control=glmerControl(optCtrl=list(maxfun=2e4)))
	}
	if (length(tmpfit.M3@optinfo$conv$lme4)!=0){
	  cat("\tre-fitting model M3", "\n")
		tmp.ss.M3 <- getME(tmpfit.M3,c("theta","fixef"))
		tmpfit.M3 <- update(tmpfit.M3, start=tmp.ss.M3, control=glmerControl(optCtrl=list(maxfun=2e4)))
	}
	
  # if any still not converged then mark this sim as not converged
	if (any(!is.null(tmpfit.M0@optinfo$conv$lme4$messages), !is.null(tmpfit.M1@optinfo$conv$lme4$messages), 
	        !is.null(tmpfit.M2@optinfo$conv$lme4$messages), !is.null(tmpfit.M3@optinfo$conv$lme4$messages) ) ) {
		simout.ulgr.m3$converge[k]=FALSE		
	}

  # collect outputs
  simout.ulgr.m3[k, 1:4] <- summary(tmpfit.M1)$coef["nbypass", ]
  simout.ulgr.m3[k, 5] <- AIC(tmpfit.M1)
  simout.ulgr.m3[k, 6:9] <- summary(tmpfit.M2)$coef["nbypass", ]
  simout.ulgr.m3[k, 10] <- AIC(tmpfit.M2)
  simout.ulgr.m3[k, 11] <- AIC(tmpfit.M3)
  simout.ulgr.m3[k, 12] <- AIC(tmpfit.M0)

  write.csv(simout.ulgr.m3, file="output/chinook_ulgr_dataM3_simulation_outputs.csv", quote=FALSE, row.names=FALSE)
  if (k==nsim) cat("Done!", "\n")
}



### ----- Summarize outputs -------------------------

## remove any non-converged runs 
simout.ulgr.m3f <- simout.ulgr.m3
simout.ulgr.m3 <- simout.ulgr.m3[simout.ulgr.m3$converge==TRUE & !is.na(simout.ulgr.m3$est.nbyp.M1), ]
# check if there are at least 1000 remaining sims first
dim(simout.ulgr.m3)

## reduce to first 1000 of remaining and run additional if necessary
nsim <- 1000
simout.ulgr.m3 <- simout.ulgr.m3[1:nsim, ]

## create outputs for Table 3
p.n.ulg.m3.m1 <- round(100*sum(simout.ulgr.m3$est.nbyp.M1 < 0 )/nsim, 1)
p.ns.ulg.m3.m1 <- round(100*sum(simout.ulgr.m3$est.nbyp.M1 < 0 & simout.ulgr.m3$pval.nbyp.M1 <= 0.05 )/nsim, 1)
p.na.ulg.m3.m1 <- round(100*sum(simout.ulgr.m3$est.nbyp.M1 < 0 & simout.ulgr.m3$aic.M1 < simout.ulgr.m3$aic.M0)/nsim, 1)

p.n.ulg.m3.m2 <- round(100*sum(simout.ulgr.m3$est.nbyp.M2 < 0)/nsim, 1)
p.ns.ulg.m3.m2 <- round(100*sum(simout.ulgr.m3$est.nbyp.M2 < 0 & simout.ulgr.m3$pval.nbyp.M2 <= 0.05 )/nsim, 1)
p.na.ulg.m3.m2 <- round(100*sum(simout.ulgr.m3$est.nbyp.M2 < 0 & simout.ulgr.m3$aic.M2 < simout.ulgr.m3$aic.M3)/nsim, 1)

outcome.labs <- c("nbyp.neg", "nbyp.neg & p<0.05", "nbyp.neg & lower AIC")
expected.labs <- c(50.0, 2.5, 7.8)
tab3.ulgr.m3 <- data.frame(data.model=rep("M3", 6), fit.model=rep(c("M1","M2"),rep(3,2)), 
                    outcome=rep(outcome.labs,2), expected=rep(expected.labs,2),
                    chinook.ULGD=c(p.n.ulg.m3.m1, p.ns.ulg.m3.m1, p.na.ulg.m3.m1, p.n.ulg.m3.m2, p.ns.ulg.m3.m2, p.na.ulg.m3.m2  ) )

tab3.ulgr <- rbind(tab3.ulgr.m0, tab3.ulgr.m3)
tab3.ulgr

## summaries of parameter estimates for Figure A1
summary(simout.ulgr.m3$est.nbyp.M1)
summary(simout.ulgr.m3$est.nbyp.M2)

quantile(simout.ulgr.m3$est.nbyp.M1, probs=c(0.025, 0.975) )
quantile(simout.ulgr.m3$est.nbyp.M2, probs=c(0.025, 0.975) )




####---------------------------------------------------------
##     Tagged at Lower Granite
####---------------------------------------------------------


###-----------------------------
##     Generate Data from M0
###-----------------------------


## ---- Set up covariance matrices and parameter means------
rcor.lgr.m0 <-  as.data.frame(VarCorr(modfits.lgr[[1]]))
G1.lgr.m0 <- matrix(0, 22, 22)
diag(G1.lgr.m0)[seq(1,21,by=2)] <- rcor.lgr.m0$vcov[1]
diag(G1.lgr.m0)[seq(2,22,by=2)] <- rcor.lgr.m0$vcov[2]
for (j in seq(1,21,by=2)){
	G1.lgr.m0[j, j+1] <- G1.lgr.m0[j+1, j] <- rcor.lgr.m0$vcov[3]
}
G.lgr.m0 <- as(G1.lgr.m0, Class="Matrix")  
X.lgr.m0 <- as(getME(modfits.lgr[[1]], "X"), Class="Matrix")						
Z.lgr.m0 <- as(getME(modfits.lgr[[1]], "Z"), Class="Matrix")						 
VC.lgr.m0 <- fullCov(modfits.lgr[[1]], G.lgr.m0)
bmu.lgr.m0 <- fixef(modfits.lgr[[1]])
rand.lgr.m0 <- as.matrix(ranef(modfits.lgr[[1]])$fyear)
rmu.lgr.m0 <- interleave(rand.lgr.m0[,1], rand.lgr.m0[,2])



### --- Run simulations ------------------------

formula.lgr.M0 <- formula(paste("adult_return ~ ", covsets.lgr[[bestmod.cov.lgr]] ))

nsim <- 1200

simout.lgr.m0 <- data.frame(matrix(NA, nsim, 7 ))
names(simout.lgr.m0) <- c("est.nbyp.M1", "SE.nbyp.M1", "zval.nbyp.M1", "pval.nbyp.M1","aic.M1", "aic.M0", "converge")
simout.lgr.m0$converge=TRUE

tmpdat.lgr <- cdat.lgr

for (k in 1:nsim) {
  set.seed(k*20) 
  cat("sim ", k, " of", nsim, "\n")
  tmp.pars <- rmvnorm(n=1, mean=c(bmu.lgr.m0, rmu.lgr.m0), sigma=VC.lgr.m0 )
  tmp.pred.eta <- X.lgr.m0%*%as(as.matrix(tmp.pars[1:3]), Class="Matrix") + Z.lgr.m0%*%as(as.matrix(tmp.pars[4:25]), Class="Matrix")
  tmp.pred <- plogis(as.vector(as(tmp.pred.eta, Class="matrix")))
  # generate new adult return data
  tmp.y <- rbinom(n=length(tmp.pred.eta), size=1, prob=tmp.pred)
  tmpdat.lgr$adult_return <- tmp.y

  # fit models
  cat("\tfitting model M0", "\n")
  tmpfit.M0 <- glmer(formula.lgr.M0, family=binomial, data=tmpdat.lgr)
  cat("\tfitting model M1", "\n")
  tmpfit.M1 <- glmer(formulas.lgr[[4]], family=binomial, data=tmpdat.lgr)

  # re-fit any that did not converge
  if (length(tmpfit.M0@optinfo$conv$lme4)!=0){
    cat("\tre-fitting model M0", "\n")
    tmp.ss.M0 <- getME(tmpfit.M0,c("theta","fixef"))
    tmpfit.M0 <- update(tmpfit.M0, start=tmp.ss.M0,control=glmerControl(optCtrl=list(maxfun=2e4)))
  }
  if (length(tmpfit.M1@optinfo$conv$lme4)!=0){
    cat("\tre-fitting model M1", "\n")
    tmp.ss.M1 <- getME(tmpfit.M1,c("theta","fixef"))
    tmpfit.M1 <- update(tmpfit.M1, start=tmp.ss.M1,control=glmerControl(optCtrl=list(maxfun=2e4)))
  }
	
  # if any still not converged then mark this sim as not converged
  if (any(!is.null(tmpfit.M0@optinfo$conv$lme4$messages), !is.null(tmpfit.M1@optinfo$conv$lme4$messages)) ) {
    simout.lgr.m0$converge[k]=FALSE		
  }

  # collect outputs
  simout.lgr.m0[k, 1:4] <- summary(tmpfit.M1)$coef["nbypass", ]
  simout.lgr.m0[k, 5] <- AIC(tmpfit.M1)
  simout.lgr.m0[k, 6] <- AIC(tmpfit.M0)

  write.csv(simout.lgr.m0, file="output/chinook_lgr_dataM0_simulation_outputs.csv", quote=FALSE, row.names=FALSE)
  if (k==nsim) cat("Done!", "\n")
}



### ----- Summarize outputs -------------------------

## remove any non-converged runs 
simout.lgr.m0f <- simout.lgr.m0
simout.lgr.m0 <- simout.lgr.m0[simout.lgr.m0$converge==TRUE & !is.na(simout.lgr.m0$est.nbyp.M1), ]
# check if there are at least 1000 remaining sims first
dim(simout.lgr.m0)

## reduce to first 1000 of remaining and run additional if necessary
nsim <- 1000
simout.lgr.m0 <- simout.lgr.m0[1:nsim, ]

## create outputs for Table 3

# nsim <- 332
p.n.lg.m0.m1 <- round(100*sum(simout.lgr.m0$est.nbyp.M1 < 0 )/nsim, 1)
p.ns.lg.m0.m1 <- round(100*sum(simout.lgr.m0$est.nbyp.M1 < 0 & simout.lgr.m0$pval.nbyp.M1 <= 0.05 )/nsim, 1)
p.na.lg.m0.m1 <- round(100*sum(simout.lgr.m0$est.nbyp.M1 < 0 & simout.lgr.m0$aic.M1 < simout.lgr.m0$aic.M0)/nsim, 1)

outcome.labs <- c("nbyp.neg", "nbyp.neg & p<0.05", "nbyp.neg & lower AIC")
expected.labs <- c(50.0, 2.5, 7.8)
tab3.lgr.m0 <- data.frame(data.model=rep("M0", 3), fit.model=rep("M1",3), 
                    outcome=outcome.labs, expected=expected.labs,
                    chinook.LGD=c(p.n.lg.m0.m1,p.ns.lg.m0.m1,p.na.lg.m0.m1) )



###----------------------------------
##      Generate Data from M3
###----------------------------------


## ---- Set up covariance matrices and parameter means------
rcor.lgr.m3 <-  as.data.frame(VarCorr(modfits.lgr[[2]]))
G1.lgr.m3 <- matrix(0, 22, 22)
diag(G1.lgr.m3)[seq(1,21,by=2)] <- rcor.lgr.m3$vcov[1]
diag(G1.lgr.m3)[seq(2,22,by=2)] <- rcor.lgr.m3$vcov[2]
for (j in seq(1,21,by=2)){
	G1.lgr.m3[j, j+1] <- G1.lgr.m3[j+1, j] <- rcor.lgr.m3$vcov[3]
}
G.lgr.m3 <- as(G1.lgr.m3, Class="Matrix")  
X.lgr.m3 <- as(getME(modfits.lgr[[2]], "X"), Class="Matrix")						
Z.lgr.m3 <- as(getME(modfits.lgr[[2]], "Z"), Class="Matrix")						 
VC.lgr.m3 <- fullCov(modfits.lgr[[2]], G.lgr.m3)
bmu.lgr.m3 <- fixef(modfits.lgr[[2]])
rand.lgr.m3 <- as.matrix(ranef(modfits.lgr[[2]])$fyear)
rmu.lgr.m3 <- interleave(rand.lgr.m3[,1], rand.lgr.m3[,2])


### --- Run simulations ------------------------

formula.lgr.M0 <- formula(paste("adult_return ~ ", covsets.lgr[[bestmod.cov.lgr]] ))

nsim <- 1200  

simout.lgr.m3 <- data.frame(matrix(NA, nsim, 13 ))
names(simout.lgr.m3) <- c("est.nbyp.M1", "SE.nbyp.M1", "zval.nbyp.M1", "pval.nbyp.M1","aic.M1","est.nbyp.M2", "SE.nbyp.M2", "zval.nbyp.M2", 
"pval.nbyp.M2","aic.M2","aic.M3", "aic.M0",  "converge")
simout.lgr.m3$converge=TRUE

tmpdat.lgr <- cdat.lgr

for (k in 1:nsim) {
	set.seed(k*18) 
  cat("sim ", k, " of", nsim, "\n")
  tmp.pars <- rmvnorm(n=1, mean=c(bmu.lgr.m3, rmu.lgr.m3), sigma=VC.lgr.m3 )
  tmp.pred.eta <- X.lgr.m3%*%as(as.matrix(tmp.pars[1:4]), Class="Matrix") + Z.lgr.m3%*%as(as.matrix(tmp.pars[5:26]), Class="Matrix")
  tmp.pred <- plogis(as.vector(as(tmp.pred.eta, Class="matrix")))
  # generate new adult return data
  tmp.y <- rbinom(n=length(tmp.pred.eta), size=1, prob=tmp.pred)
  tmpdat.lgr$adult_return <- tmp.y

  # fit models
  cat("\tfitting model M0", "\n")
  tmpfit.M0 <- glmer(formula.lgr.M0, family=binomial, data=tmpdat.lgr)
  cat("\tfitting model M1", "\n")
  tmpfit.M1 <- glmer(formulas.lgr[[4]], family=binomial, data=tmpdat.lgr)
  cat("\tfitting model M2", "\n")
  tmpfit.M2 <- glmer(formulas.lgr[[7]], family=binomial, data=tmpdat.lgr)
  cat("\tfitting model M3", "\n")
  tmpfit.M3 <- glmer(formulas.lgr[[1]], family=binomial, data=tmpdat.lgr)

  # re-fit any that did not converge
	if (length(tmpfit.M0@optinfo$conv$lme4)!=0){
	  cat("\tre-fitting model M0", "\n")
		tmp.ss.M0 <- getME(tmpfit.M0,c("theta","fixef"))
		tmpfit.M0 <- update(tmpfit.M0, start=tmp.ss.M0,control=glmerControl(optCtrl=list(maxfun=2e4)))
	}
	if (length(tmpfit.M1@optinfo$conv$lme4)!=0){
	  cat("\tre-fitting model M1", "\n")
		tmp.ss.M1 <- getME(tmpfit.M1,c("theta","fixef"))
		tmpfit.M1 <- update(tmpfit.M1, start=tmp.ss.M1,control=glmerControl(optCtrl=list(maxfun=2e4)))
	}
	if (length(tmpfit.M2@optinfo$conv$lme4)!=0){
	  cat("\tre-fitting model M2", "\n")
		tmp.ss.M2 <- getME(tmpfit.M2,c("theta","fixef"))
		tmpfit.M2 <- update(tmpfit.M2, start=tmp.ss.M2, control=glmerControl(optCtrl=list(maxfun=2e4)))
	}
	if (length(tmpfit.M3@optinfo$conv$lme4)!=0){
	  cat("\tre-fitting model M3", "\n")
		tmp.ss.M3 <- getME(tmpfit.M3,c("theta","fixef"))
		tmpfit.M3 <- update(tmpfit.M3, start=tmp.ss.M3, control=glmerControl(optCtrl=list(maxfun=2e4)))
	}
	
  # if any still not converged then mark this sim as not converged
	if (any(!is.null(tmpfit.M0@optinfo$conv$lme4$messages), !is.null(tmpfit.M1@optinfo$conv$lme4$messages), 
	        !is.null(tmpfit.M2@optinfo$conv$lme4$messages), !is.null(tmpfit.M3@optinfo$conv$lme4$messages) ) ) {
		simout.lgr.m3$converge[k]=FALSE		
	}

  # collect outputs
  simout.lgr.m3[k, 1:4] <- summary(tmpfit.M1)$coef["nbypass", ]
  simout.lgr.m3[k, 5] <- AIC(tmpfit.M1)
  simout.lgr.m3[k, 6:9] <- summary(tmpfit.M2)$coef["nbypass", ]
  simout.lgr.m3[k, 10] <- AIC(tmpfit.M2)
  simout.lgr.m3[k, 11] <- AIC(tmpfit.M3)
  simout.lgr.m3[k, 12] <- AIC(tmpfit.M0)

  write.csv(simout.lgr.m3, file="output/chinook_lgr_dataM3_simulation_outputs.csv", quote=FALSE, row.names=FALSE)
  if (k==nsim) cat("Done!", "\n")
}



### ----- Summarize outputs -------------------------

## remove any non-converged runs 
simout.lgr.m3f <- simout.lgr.m3
simout.lgr.m3 <- simout.lgr.m3[simout.lgr.m3$converge==TRUE & !is.na(simout.lgr.m3$est.nbyp.M1), ]
# check if there are at least 1000 remaining sims first
dim(simout.lgr.m3)

## reduce to first 1000 of remaining and run additional if necessary
nsim <- 1000
simout.lgr.m3 <- simout.lgr.m3[1:nsim, ]

## create outputs for Table 3
p.n.lg.m3.m1 <- round(100*sum(simout.lgr.m3$est.nbyp.M1 < 0 )/nsim, 1)
p.ns.lg.m3.m1 <- round(100*sum(simout.lgr.m3$est.nbyp.M1 < 0 & simout.lgr.m3$pval.nbyp.M1 <= 0.05 )/nsim, 1)
p.na.lg.m3.m1 <- round(100*sum(simout.lgr.m3$est.nbyp.M1 < 0 & simout.lgr.m3$aic.M1 < simout.lgr.m3$aic.M0)/nsim, 1)

p.n.lg.m3.m2 <- round(100*sum(simout.lgr.m3$est.nbyp.M2 < 0)/nsim, 1)
p.ns.lg.m3.m2 <- round(100*sum(simout.lgr.m3$est.nbyp.M2 < 0 & simout.lgr.m3$pval.nbyp.M2 <= 0.05 )/nsim, 1)
p.na.lg.m3.m2 <- round(100*sum(simout.lgr.m3$est.nbyp.M2 < 0 & simout.lgr.m3$aic.M2 < simout.lgr.m3$aic.M3)/nsim, 1)

outcome.labs <- c("nbyp.neg", "nbyp.neg & p<0.05", "nbyp.neg & lower AIC")
expected.labs <- c(50.0, 2.5, 7.8)
tab3.lgr.m3 <- data.frame(data.model=rep("M3", 6), fit.model=rep(c("M1","M2"),rep(3,2)), 
                    outcome=rep(outcome.labs,2), expected=rep(expected.labs,2),
                    chinook.LGD=c(p.n.lg.m3.m1,p.ns.lg.m3.m1,p.na.lg.m3.m1, p.n.lg.m3.m2,p.ns.lg.m3.m2,p.na.lg.m3.m2  ) )

tab3.lgr <- rbind(tab3.lgr.m0, tab3.lgr.m3)
tab3.lgr

## summaries of parameter estimates for Figure A1
summary(simout.lgr.m3$est.nbyp.M1)
summary(simout.lgr.m3$est.nbyp.M2)

quantile(simout.lgr.m3$est.nbyp.M1, probs=c(0.025, 0.975) )
quantile(simout.lgr.m3$est.nbyp.M2, probs=c(0.025, 0.975) )




### Create combined output for Table 3

tab3 <- merge(tab3.ulgr, tab3.lgr, sort=FALSE)

write.csv(tab3, file="output/chinook_simulations_table3.csv", quote=F, row.names=F)











































