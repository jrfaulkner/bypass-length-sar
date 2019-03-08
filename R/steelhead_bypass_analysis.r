## Steelhead Bypass Probability Analysis

## For model fitting in support of:
## Faulkner, JR, BL Bellerud, DL Widener, and RW Zabel. 2019. Associations among fish length, dam passage history, and survival to adulthood in two at-risk species of Pacific salmon.


library(lme4)



###-------------------------------------------------------------------------------------
##                 Data
###-------------------------------------------------------------------------------------

### --- Read data files
sdat <- read.csv("../../data/steelhead_bypass_data.csv", stringsAsFactors=F)


## --- Make new variables
## make year a factor
sdat$fyear <- as.character(sdat$year)






###----------------------------------------------------------------------------------------
##            Model Fitting
###----------------------------------------------------------------------------------------


### ----  Functions ---------------------------------

# function to fit a set of models and output to a list
fitModels <- function(modset, respname, data, subset){
  nm <- length(modset)
  convec <- numeric(nm)
  outlist <- list()
  for (kk in 1:nm) {
    cat("fitting model: ", kk, " of ", nm, "\n" )
    tmp.form <- formula(paste(respname, "~ ", modset[kk] ))
    tmp.fit <- glmer(tmp.form, family=binomial, data=data, subset=subset, control=glmerControl(optCtrl=list(maxfun=2e5)) )
    if (length(tmp.fit@optinfo$conv$lme4)==0) convec[kk] <- 0
    if (length(tmp.fit@optinfo$conv$lme4)!=0) convec[kk] <- tmp.fit@optinfo$conv$lme4$code
    outlist[[kk]] <- tmp.fit 
    # if did not converge, try to refit from new start
    if (convec[kk]!=0) {
      cat("re-fitting model: ", kk, " of ", nm, "\n" )
      tmp.ss <- getME(tmp.fit,c("theta","fixef"))
      tmp.fit2 <- update(tmp.fit,start=tmp.ss,control=glmerControl(optCtrl=list(maxfun=2e5)) )
      if (length(tmp.fit2@optinfo$conv$lme4)==0) convec[kk] <- 0
      if (length(tmp.fit2@optinfo$conv$lme4)!=0) convec[kk] <- tmp.fit2@optinfo$conv$lme4$code
      outlist[[kk]] <- tmp.fit2 
    }
    # if still did not converge, try to refit with different optimizer
    if (convec[kk]!=0) {
      cat("re-fitting model: ", kk, " of ", nm, "\n" )
      tmp.ss2 <- getME(tmp.fit2,c("theta","fixef"))
      tmp.fit3 <- update(tmp.fit2,start=tmp.ss2,control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
      if (length(tmp.fit3@optinfo$conv$lme4)==0) convec[kk] <- 0
      if (length(tmp.fit3@optinfo$conv$lme4)!=0) convec[kk] <- tmp.fit3@optinfo$conv$lme4$code
      outlist[[kk]] <- tmp.fit3 
    }
    if (convec[kk]!=0) cat("2nd re-fit of model ", kk, " failed!", "\n" )
  }
 outlist$converge <- convec
 return(outlist)
}



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
# R = release site 
# G = day at Lower Granite 
# G2 = G squared
# W = week at Lower Granite
# D = day at Bonneville or Trawl
# D2 = D squared
# L = length


### Possible covariate-only models
# notes: all models have Y and 
#   any G, W, or D variable cannot occur together (only one set of time variables)
# All models are constructed separately by reartype and Dam

# 1) Y
# 2) Y + R
# 3) Y + G
# 4) Y + G + G2
# 5) Y + R + G
# 6) Y + R + G + G2
# 7) Y + G + Y*G
# 8) Y + G + G2 + Y*G
# 9) Y + G + G2 + Y*G + Y*G2
# 10) Y + R + G + Y*G
# 11) Y + R + G + G2 + Y*G
# 12) Y + R + G + G2 + Y*G + Y*G2
# 13) Y + D
# 14) Y + D + D2
# 15) Y + R + D
# 16) Y + R + D + D2
# 17) Y + D + Y*D
# 18) Y + D + D2 + Y*D
# 19) Y + D + D2 + Y*D + Y*D2
# 20) Y + R + D + Y*D
# 21) Y + R + D + D2 + Y*D
# 22) Y + R + D + D2 + Y*D + Y*D2
# 23) Y + W
# 24) Y + R + W
# 25) Y + W + Y*W
# 26) Y + R + W + Y*W




covsets <- c("(1 | fyear)",
             "rel_site2 + (1 | fyear)",
             "zday_lgr + (1 | fyear)",
             "zday_lgr + zday_lgr_sq + (1 | fyear)",
             "rel_site2 + zday_lgr + (1 | fyear)",
             "rel_site2 + zday_lgr + zday_lgr_sq + (1 | fyear)", 
             "zday_lgr + (zday_lgr | fyear)",
             "zday_lgr + zday_lgr_sq + (zday_lgr | fyear)",
             "zday_lgr + zday_lgr_sq + (zday_lgr + zday_lgr_sq | fyear)", 
             "rel_site2 + zday_lgr + (zday_lgr | fyear)",
             "rel_site2 + zday_lgr + zday_lgr_sq + (zday_lgr | fyear)",
             "rel_site2 + zday_lgr + zday_lgr_sq + (zday_lgr + zday_lgr_sq | fyear)",  
             "zlast_day + (1 | fyear)",
             "zlast_day + zlast_day_sq + (1 | fyear)",
             "rel_site2 + zlast_day + (1 | fyear)",
             "rel_site2 + zlast_day + zlast_day_sq + (1 | fyear)", 
             "zlast_day + (zlast_day | fyear)",
             "zlast_day + zlast_day_sq + (zlast_day | fyear)",
             "zlast_day + zlast_day_sq + (zlast_day + zlast_day_sq | fyear)", 
             "rel_site2 + zlast_day + (zlast_day | fyear)",
             "rel_site2 + zlast_day + zlast_day_sq + (zlast_day | fyear)",
             "rel_site2 + zlast_day + zlast_day_sq + (zlast_day + zlast_day_sq | fyear)", 
             "rel_week_lgr + (1 | fyear)",
             "rel_site2 + rel_week_lgr + (1 | fyear)",
             "rel_week_lgr + (rel_week_lgr | fyear)",
             "rel_site2 + rel_week_lgr + (rel_week_lgr | fyear)" )

cov.modnums <- 1:length(covsets)

testsets <- "zlength"

respsets <- c("bypass_LGR", "bypass_LGS", "bypass_LMN", "bypass_IHR", "bypass_MCN", "bypass_JDA", "bypass_BON", "bypass_BON")

damsets <- c("LGR", "LGS", "LMN", "IHR", "MCN", "JDA", "BON", "BCC")

subsets <- list(sdat$rel_site2!="LGR", sdat$last_site=="twx", (sdat$bypass_BON==1 | sdat$pass_BCC==1) & sdat$fyear %in% as.character(2006:2014))


nmods.cov <- length(covsets)



###-----------------------------------------------------------------------
##               Wild Fish (unclipped)
###-----------------------------------------------------------------------

###--------------------------------------------------
##        Lower Granite Dam 
###--------------------------------------------------

idam <- which(damsets=="LGR")
# -------   Fit covariate-only models -------------------------------------

modfits.cov.lgr.w <- fitModels(covsets, respsets[idam], sdat, subsets[[1]] & sdat$clipped=="b. no")

# check if any models did not converge (any non-zero)
modfits.cov.lgr.w$converge
# all are fine

# Calculate AIC
aicsum.cov.lgr.w <- getAICsummary(modfits.cov.lgr.w, covsets)
aicsum.cov.lgr.w

# find best model
bestmod.cov.lgr.w <- aicsum.cov.lgr.w$Model[1]

# -------   Fit model with length -------------------------------------

tmpset <- paste("zlength +", covsets[bestmod.cov.lgr.w])
modfit.len.lgr.w <- fitModels(tmpset, respsets[idam], sdat, subsets[[1]] & sdat$clipped=="b. no")
modfit.len.lgr.w$converge

modsum.len.lgr.w <- summary(modfit.len.lgr.w[[1]])
modsum.len.lgr.w






###--------------------------------------------------
##        Little Goose Dam 
###--------------------------------------------------

idam <- which(damsets=="LGS")
# -------   Fit covariate-only models -------------------------------------
modfits.cov.lgs.w <- fitModels(covsets, respsets[idam], sdat, sdat$clipped=="b. no")

# check if any models did not converge (any non-zero)
modfits.cov.lgs.w$converge
# all are fine

# Calculate AIC
aicsum.cov.lgs.w <- getAICsummary(modfits.cov.lgs.w, covsets)
aicsum.cov.lgs.w

# find best model
bestmod.cov.lgs.w <- aicsum.cov.lgs.w$Model[1]

# -------   Fit model with length -------------------------------------
tmpset <- paste("zlength +", covsets[bestmod.cov.lgs.w])
modfit.len.lgs.w <- fitModels(tmpset, respsets[idam], sdat, sdat$clipped=="b. no")
modfit.len.lgs.w$converge
modsum.len.lgs.w <- summary(modfit.len.lgs.w[[1]])
modsum.len.lgs.w




###--------------------------------------------------
##        Lower Monumental Dam 
###--------------------------------------------------

idam <- which(damsets=="LMN")
# -------   Fit covariate-only models -------------------------------------
modfits.cov.lmn.w <- fitModels(covsets, respsets[idam], sdat, sdat$clipped=="b. no")

# check if any models did not converge (any non-zero)
modfits.cov.lmn.w$converge
# all are fine

# Calculate AIC
aicsum.cov.lmn.w <- getAICsummary(modfits.cov.lmn.w, covsets)
aicsum.cov.lmn.w

# find best model
bestmod.cov.lmn.w <- aicsum.cov.lmn.w$Model[1]

# -------   Fit model with length -------------------------------------

tmpset <- paste("zlength +", covsets[bestmod.cov.lmn.w])
modfit.len.lmn.w <- fitModels(tmpset, respsets[idam], sdat, sdat$clipped=="b. no")
modfit.len.lmn.w$converge
modsum.len.lmn.w <- summary(modfit.len.lmn.w[[1]])
modsum.len.lmn.w




###--------------------------------------------------
##        Ice Harbor Dam 
###--------------------------------------------------

idam <- which(damsets=="IHR")
# -------   Fit covariate-only models -------------------------------------
modfits.cov.ihr.w <- fitModels(covsets, respsets[idam], sdat, sdat$clipped=="b. no" & sdat$fyear %in% as.character(2005:2014))

# check if any models did not converge (any non-zero)
modfits.cov.ihr.w$converge 
# all are fine 

# Calculate AIC
aicsum.cov.ihr.w <- getAICsummary(modfits.cov.ihr.w, covsets)
aicsum.cov.ihr.w

# find best model
bestmod.cov.ihr.w <- aicsum.cov.ihr.w$Model[1]

# -------   Fit model with length -------------------------------------

tmpset <- paste("zlength +", covsets[bestmod.cov.ihr.w])
modfit.len.ihr.w <- fitModels(tmpset, respsets[idam], sdat, sdat$clipped=="b. no" & sdat$fyear %in% as.character(2005:2014))
modfit.len.ihr.w$converge
modsum.len.ihr.w <- summary(modfit.len.ihr.w[[1]])
modsum.len.ihr.w



###--------------------------------------------------
##        McNary Dam 
###--------------------------------------------------

idam <- which(damsets=="MCN")
# -------   Fit covariate-only models -------------------------------------
modfits.cov.mcn.w <- fitModels(covsets, respsets[idam], sdat, sdat$clipped=="b. no")

# check if any models did not converge (any non-zero)
modfits.cov.mcn.w$converge 
# all are fine

# Calculate AIC
aicsum.cov.mcn.w <- getAICsummary(modfits.cov.mcn.w, covsets)
aicsum.cov.mcn.w

# find best model
bestmod.cov.mcn.w <- aicsum.cov.mcn.w$Model[1]

# -------   Fit model with length -------------------------------------

tmpset <- paste("zlength +", covsets[bestmod.cov.mcn.w])
modfit.len.mcn.w <- fitModels(tmpset, respsets[idam], sdat, sdat$clipped=="b. no")
modfit.len.mcn.w$converge
modsum.len.mcn.w <- summary(modfit.len.mcn.w[[1]])
modsum.len.mcn.w


###--------------------------------------------------
##        John Day Dam 
###--------------------------------------------------

idam <- which(damsets=="JDA")
# -------   Fit covariate-only models -------------------------------------
modfits.cov.jda.w <- fitModels(covsets, respsets[idam], sdat, sdat$clipped=="b. no")

# check if any models did not converge (any non-zero)
modfits.cov.jda.w$converge
# all are fine

# Calculate AIC
aicsum.cov.jda.w <- getAICsummary(modfits.cov.jda.w, covsets)
aicsum.cov.jda.w

# find best model
bestmod.cov.jda.w <- aicsum.cov.jda.w$Model[1]

# -------   Fit model with length -------------------------------------

tmpset <- paste("zlength +", covsets[bestmod.cov.jda.w])
modfit.len.jda.w <- fitModels(tmpset, respsets[idam], sdat, sdat$clipped=="b. no")
modfit.len.jda.w$converge
modsum.len.jda.w <- summary(modfit.len.jda.w[[1]])
modsum.len.jda.w



###--------------------------------------------------
##        Bonneville Dam 
###--------------------------------------------------

idam <- which(damsets=="BON")
# -------   Fit covariate-only models -------------------------------------
modfits.cov.bon.w <- fitModels(covsets, respsets[idam], sdat, subsets[[2]] & sdat$clipped=="b. no")

# check if any models did not converge (any non-zero)
modfits.cov.bon.w$converge
# all are fine

# Calculate AIC
aicsum.cov.bon.w <- getAICsummary(modfits.cov.bon.w, covsets)
aicsum.cov.bon.w

# find best model
bestmod.cov.bon.w <- aicsum.cov.bon.w$Model[1]

# -------   Fit model with length -------------------------------------

tmpset <- paste("zlength +", covsets[bestmod.cov.bon.w])
modfit.len.bon.w <- fitModels(tmpset, respsets[idam], sdat, subsets[[2]] & sdat$clipped=="b. no")
modfit.len.bon.w$converge
modsum.len.bon.w <- summary(modfit.len.bon.w[[1]])
modsum.len.bon.w




###--------------------------------------------------
##        Bonneville Dam Corner Collector
###--------------------------------------------------


idam <- which(damsets=="BCC")
# -------   Fit covariate-only models -------------------------------------
modfits.cov.bcc.w <- fitModels(covsets, respsets[idam], sdat, subsets[[3]] & sdat$clipped=="b. no")

# check if any models did not converge (any non-zero)
modfits.cov.bcc.w$converge 
# note that last model did not converge

modfits.cov.bcc.w <- modfits.cov.bcc.w[-26]

# Calculate AIC
aicsum.cov.bcc.w <- getAICsummary(modfits.cov.bcc.w, covsets[-26])
aicsum.cov.bcc.w

# find best model 
# NOTE: the best model does not converge when length is added. Therefore use #2 model
bestmod.cov.bcc.w <- aicsum.cov.bcc.w$Model[1]

# -------   Fit model with length -------------------------------------

tmpset <- paste("zlength +", covsets[bestmod.cov.bcc.w])
modfit.len.bcc.w <- fitModels(tmpset, respsets[idam], sdat, subsets[[3]] & sdat$clipped=="b. no")
modfit.len.bcc.w$converge
modsum.len.bcc.w <- summary(modfit.len.bcc.w[[1]])
modsum.len.bcc.w





## -----------------------------------------
##     Create Summary for Table 1
## -----------------------------------------

# species, reartype, dam, n, estimate, 95%CI, dAIC, p-value

dnames <- c("LGD","LGSD","LMD","IHD","MCD","JDD","BVD","BVCC")

tab1.w <- data.frame(Species=rep("Steelhead",8), Type=rep("W", 8), Dam=dnames, 
  n = c(length(modsum.len.lgr.w$resid),
    length(modsum.len.lgs.w$resid),
    length(modsum.len.lmn.w$resid),
    length(modsum.len.ihr.w$resid),
    length(modsum.len.mcn.w$resid),
    length(modsum.len.jda.w$resid),
    length(modsum.len.bon.w$resid),
    length(modsum.len.bcc.w$resid) ),
  round(rbind(extractCI(modsum.len.lgr.w, "zlength", scale="log.odds"),
    extractCI(modsum.len.lgs.w, "zlength", scale="log.odds"),
    extractCI(modsum.len.lmn.w, "zlength", scale="log.odds"),
    extractCI(modsum.len.ihr.w, "zlength", scale="log.odds"),
    extractCI(modsum.len.mcn.w, "zlength", scale="log.odds"),
    extractCI(modsum.len.jda.w, "zlength", scale="log.odds"),
    extractCI(modsum.len.bon.w, "zlength", scale="log.odds"),
    extractCI(modsum.len.bcc.w, "zlength", scale="log.odds") ),3 ),
  dAIC=c(round(modsum.len.lgr.w$AIC[1] - aicsum.cov.lgr.w$AIC[1], 1),
    round(modsum.len.lgs.w$AIC[1] - aicsum.cov.lgs.w$AIC[1], 1),
    round(modsum.len.lmn.w$AIC[1] - aicsum.cov.lmn.w$AIC[1], 1),
    round(modsum.len.ihr.w$AIC[1] - aicsum.cov.ihr.w$AIC[1], 1),
    round(modsum.len.mcn.w$AIC[1] - aicsum.cov.mcn.w$AIC[1], 1),
    round(modsum.len.jda.w$AIC[1] - aicsum.cov.jda.w$AIC[1], 1),
    round(modsum.len.bon.w$AIC[1] - aicsum.cov.bon.w$AIC[1], 1),
    round(modsum.len.bcc.w$AIC[1] - aicsum.cov.bcc.w$AIC[1], 1) ) ,
  pvalue=c(round(modsum.len.lgr.w$coef["zlength", 4], 4),
    round(modsum.len.lgs.w$coef["zlength", 4], 4),
    round(modsum.len.lmn.w$coef["zlength", 4], 4),
    round(modsum.len.ihr.w$coef["zlength", 4], 4),
    round(modsum.len.mcn.w$coef["zlength", 4], 4),
    round(modsum.len.jda.w$coef["zlength", 4], 4),
    round(modsum.len.bon.w$coef["zlength", 4], 4),
    round(modsum.len.bcc.w$coef["zlength", 4], 4) ) 
)

names(tab1.w)[5:7] <- c("Estimate", "95%CI.L", "95%CI.U")

tab1.w



## ------------------------------------------
##     Create Summary for Table A3 
## ------------------------------------------

tabA3.w <- data.frame(Species=rep("Steelhead",8), Type=rep("W", 8), Dam=dnames, 
  modform=c(aicsum.cov.lgr.w$ModelForm[1],
    aicsum.cov.lgs.w$ModelForm[1],
    aicsum.cov.lmn.w$ModelForm[1],
    aicsum.cov.ihr.w$ModelForm[1],
    aicsum.cov.mcn.w$ModelForm[1],
    aicsum.cov.jda.w$ModelForm[1],
    aicsum.cov.bon.w$ModelForm[1],
    aicsum.cov.bcc.w$ModelForm[1])
)

tabA3.w



###----------------------------------------------------------------------------------------------------
##               Hatchery Fish (clipped)
###----------------------------------------------------------------------------------------------------

###--------------------------------------------------
##        Lower Granite Dam 
###--------------------------------------------------

idam <- which(damsets=="LGR")

# -------   Fit covariate-only models -------------------------------------
modfits.cov.lgr.h <- fitModels(covsets, respsets[idam], sdat, subsets[[1]] & sdat$clipped=="a. yes")

# check if any models did not converge (any non-zero)
modfits.cov.lgr.h$converge  
# all are fine

# Calculate AIC
aicsum.cov.lgr.h <- getAICsummary(modfits.cov.lgr.h, covsets)
aicsum.cov.lgr.h

# find best model
bestmod.cov.lgr.h <- aicsum.cov.lgr.h$Model[1]

# -------   Fit model with length -------------------------------------

tmpset <- paste("zlength + ", covsets[bestmod.cov.lgr.h])
modfit.len.lgr.h <- fitModels(tmpset, respsets[idam], sdat, subsets[[1]] & sdat$clipped=="a. yes")
modfit.len.lgr.h$converge 

modsum.len.lgr.h <- summary(modfit.len.lgr.h[[1]])
modsum.len.lgr.h



###--------------------------------------------------
##        Little Goose Dam 
###--------------------------------------------------

idam <- which(damsets=="LGS")
# -------   Fit covariate-only models -------------------------------------
modfits.cov.lgs.h <- fitModels(covsets, respsets[idam], sdat, sdat$clipped=="a. yes")

# check if any models did not converge (any non-zero)
modfits.cov.lgs.h$converge 
# all are fine

# Calculate AIC
aicsum.cov.lgs.h <- getAICsummary(modfits.cov.lgs.h, covsets)
aicsum.cov.lgs.h

# find best model
bestmod.cov.lgs.h <- aicsum.cov.lgs.h$Model[1]

# -------   Fit model with length -------------------------------------

tmpset <- paste("zlength +", covsets[bestmod.cov.lgs.h])
modfit.len.lgs.h <- fitModels(tmpset, respsets[idam], sdat, sdat$clipped=="a. yes")
modfit.len.lgs.h$converge

modsum.len.lgs.h <- summary(modfit.len.lgs.h[[1]])
modsum.len.lgs.h




###--------------------------------------------------
##        Lower Monumental Dam 
###--------------------------------------------------

idam <- which(damsets=="LMN")
# -------   Fit covariate-only models -------------------------------------
modfits.cov.lmn.h <- fitModels(covsets, respsets[idam], sdat, sdat$clipped=="a. yes")

# check if any models did not converge (any non-zero)
modfits.cov.lmn.h$converge 
# all are fine

# Calculate AIC
aicsum.cov.lmn.h <- getAICsummary(modfits.cov.lmn.h, covsets)
aicsum.cov.lmn.h

# find best model
bestmod.cov.lmn.h <- aicsum.cov.lmn.h$Model[1]

# -------   Fit model with length -------------------------------------

tmpset <- paste("zlength +", covsets[bestmod.cov.lmn.h])
modfit.len.lmn.h <- fitModels(tmpset, respsets[idam], sdat, sdat$clipped=="a. yes")
modfit.len.lmn.h$converge
modsum.len.lmn.h <- summary(modfit.len.lmn.h[[1]])
modsum.len.lmn.h




###--------------------------------------------------
##        Ice Harbor Dam 
###--------------------------------------------------

idam <- which(damsets=="IHR")
# -------   Fit covariate-only models -------------------------------------
modfits.cov.ihr.h <- fitModels(covsets, respsets[idam], sdat, sdat$clipped=="a. yes" & sdat$fyear %in% as.character(2005:2014))

# check if any models did not converge (any non-zero)
modfits.cov.ihr.h$converge
# all are fine

# Calculate AIC
aicsum.cov.ihr.h <- getAICsummary(modfits.cov.ihr.h, covsets)
aicsum.cov.ihr.h

# find best model
bestmod.cov.ihr.h <- aicsum.cov.ihr.h$Model[1]

# -------   Fit model with length -------------------------------------

tmpset <- paste("zlength +", covsets[bestmod.cov.ihr.h])
modfit.len.ihr.h <- fitModels(tmpset, respsets[idam], sdat, sdat$clipped=="a. yes" & sdat$fyear %in% as.character(2005:2014))
modfit.len.ihr.h$converge

modsum.len.ihr.h <- summary(modfit.len.ihr.h[[1]])
modsum.len.ihr.h





###--------------------------------------------------
##        McNary Dam 
###--------------------------------------------------

idam <- which(damsets=="MCN")
# -------   Fit covariate-only models -------------------------------------
modfits.cov.mcn.h <- fitModels(covsets, respsets[idam], sdat, sdat$clipped=="a. yes")

# check if any models did not converge (any non-zero)
modfits.cov.mcn.h$converge 
# all are fine

# Calculate AIC
aicsum.cov.mcn.h <- getAICsummary(modfits.cov.mcn.h, covsets)
aicsum.cov.mcn.h

# find best model
bestmod.cov.mcn.h <- aicsum.cov.mcn.h$Model[1]

# -------   Fit model with length -------------------------------------

tmpset <- paste("zlength +", covsets[bestmod.cov.mcn.h])
modfit.len.mcn.h <- fitModels(tmpset, respsets[idam], sdat, sdat$clipped=="a. yes")
modfit.len.mcn.h$converge

modsum.len.mcn.h <- summary(modfit.len.mcn.h[[1]])
modsum.len.mcn.h


###--------------------------------------------------
##        John Day Dam 
###--------------------------------------------------

idam <- which(damsets=="JDA")
# -------   Fit covariate-only models -------------------------------------
modfits.cov.jda.h <- fitModels(covsets, respsets[idam], sdat, sdat$clipped=="a. yes")

# check if any models did not converge (any non-zero)
modfits.cov.jda.h$converge 
# all are fine

# Calculate AIC
aicsum.cov.jda.h <- getAICsummary(modfits.cov.jda.h, covsets)
aicsum.cov.jda.h

# find best model
bestmod.cov.jda.h <- aicsum.cov.jda.h$Model[1]

# -------   Fit model with length -------------------------------------
tmpset <- paste("zlength +", covsets[bestmod.cov.jda.h])
modfit.len.jda.h <- fitModels(tmpset, respsets[idam], sdat, sdat$clipped=="a. yes")
modfit.len.jda.h$converge

modsum.len.jda.h <- summary(modfit.len.jda.h[[1]])
modsum.len.jda.h



###--------------------------------------------------
##        Bonneville Dam 
###--------------------------------------------------

idam <- which(damsets=="BON")
# -------   Fit covariate-only models -------------------------------------
modfits.cov.bon.h <- fitModels(covsets, respsets[idam], sdat, subsets[[2]] & sdat$clipped=="a. yes")

# check if any models did not converge (any non-zero)
modfits.cov.bon.h$converge 
# all are fine

# Calculate AIC
aicsum.cov.bon.h <- getAICsummary(modfits.cov.bon.h, covsets)
aicsum.cov.bon.h

# find best model
bestmod.cov.bon.h <- aicsum.cov.bon.h$Model[1]

# -------   Fit model with length -------------------------------------

tmpset <- paste("zlength +", covsets[bestmod.cov.bon.h])
modfit.len.bon.h <- fitModels(tmpset, respsets[idam], sdat, subsets[[2]] & sdat$clipped=="a. yes")
modfit.len.bon.h$converge

modsum.len.bon.h <- summary(modfit.len.bon.h[[1]])
modsum.len.bon.h




###--------------------------------------------------
##        Bonneville Dam Corner Collector
###--------------------------------------------------

idam <- which(damsets=="BCC")
# -------   Fit covariate-only models -------------------------------------
modfits.cov.bcc.h <- fitModels(covsets, respsets[idam], sdat, subsets[[3]] & sdat$clipped=="a. yes")

# check if any models did not converge (any non-zero)
modfits.cov.bcc.h$converge
# all are fine

# Calculate AIC
aicsum.cov.bcc.h <- getAICsummary(modfits.cov.bcc.h, covsets)
aicsum.cov.bcc.h

# find best model
bestmod.cov.bcc.h <- aicsum.cov.bcc.h$Model[1]

# -------   Fit model with length -------------------------------------
tmpset <- paste("zlength +", covsets[bestmod.cov.bcc.h])
modfit.len.bcc.h <- fitModels(tmpset, respsets[idam], sdat, subsets[[3]] & sdat$clipped=="a. yes")
modfit.len.bcc.h$converge

modsum.len.bcc.h <- summary(modfit.len.bcc.h[[1]])
modsum.len.bcc.h


## -----------------------------------------
##     Create Summary for Table 1
## -----------------------------------------

# species, reartype, dam, n, estimate, 95%CI, dAIC, p-value

dnames <- c("LGD","LGSD","LMD","IHD","MCD","JDD","BVD","BVCC")

tab1.h <- data.frame(Species=rep("Steelhead",8), Type=rep("H", 8), Dam=dnames, 
  n = c(length(modsum.len.lgr.h$resid),
    length(modsum.len.lgs.h$resid),
    length(modsum.len.lmn.h$resid),
    length(modsum.len.ihr.h$resid),
    length(modsum.len.mcn.h$resid),
    length(modsum.len.jda.h$resid),
    length(modsum.len.bon.h$resid),
    length(modsum.len.bcc.h$resid) ),
  round(rbind(extractCI(modsum.len.lgr.h, "zlength", scale="log.odds"),
    extractCI(modsum.len.lgs.h, "zlength", scale="log.odds"),
    extractCI(modsum.len.lmn.h, "zlength", scale="log.odds"),
    extractCI(modsum.len.ihr.h, "zlength", scale="log.odds"),
    extractCI(modsum.len.mcn.h, "zlength", scale="log.odds"),
    extractCI(modsum.len.jda.h, "zlength", scale="log.odds"),
    extractCI(modsum.len.bon.h, "zlength", scale="log.odds"),
    extractCI(modsum.len.bcc.h, "zlength", scale="log.odds") ),3 ),
  dAIC=c(round(modsum.len.lgr.h$AIC[1] - aicsum.cov.lgr.h$AIC[1], 1),
    round(modsum.len.lgs.h$AIC[1] - aicsum.cov.lgs.h$AIC[1], 1),
    round(modsum.len.lmn.h$AIC[1] - aicsum.cov.lmn.h$AIC[1], 1),
    round(modsum.len.ihr.h$AIC[1] - aicsum.cov.ihr.h$AIC[1], 1),
    round(modsum.len.mcn.h$AIC[1] - aicsum.cov.mcn.h$AIC[1], 1),
    round(modsum.len.jda.h$AIC[1] - aicsum.cov.jda.h$AIC[1], 1),
    round(modsum.len.bon.h$AIC[1] - aicsum.cov.bon.h$AIC[1], 1),
    round(modsum.len.bcc.h$AIC[1] - aicsum.cov.bcc.h$AIC[1], 1) ) ,
  pvalue=c(round(modsum.len.lgr.h$coef["zlength", 4], 4),
    round(modsum.len.lgs.h$coef["zlength", 4], 4),
    round(modsum.len.lmn.h$coef["zlength", 4], 4),
    round(modsum.len.ihr.h$coef["zlength", 4], 4),
    round(modsum.len.mcn.h$coef["zlength", 4], 4),
    round(modsum.len.jda.h$coef["zlength", 4], 4),
    round(modsum.len.bon.h$coef["zlength", 4], 4),
    round(modsum.len.bcc.h$coef["zlength", 4], 4) ) 
)

names(tab1.h)[5:7] <- c("Estimate", "95%CI.L", "95%CI.U")

tab1.h


tab1 <- rbind(tab1.h, tab1.w)

write.csv(tab1, file="output/steelhead_Table1_values.csv", quote=F, row.names=F)


## ------------------------------------------
##     Create Summary for Table A3 
## ------------------------------------------

tabA3.h <- data.frame(Species=rep("Steelhead",8), Type=rep("H", 8), Dam=dnames, 
  modform=c(aicsum.cov.lgr.h$ModelForm[1],
    aicsum.cov.lgs.h$ModelForm[1],
    aicsum.cov.lmn.h$ModelForm[1],
    aicsum.cov.ihr.h$ModelForm[1],
    aicsum.cov.mcn.h$ModelForm[1],
    aicsum.cov.jda.h$ModelForm[1],
    aicsum.cov.bon.h$ModelForm[1],
    aicsum.cov.bcc.h$ModelForm[1])
)

tabA3.h

tabA3 <- rbind(tabA3.h, tabA3.w)


write.csv(tabA3, file="output/steelhead_TableA3_values.csv", quote=F, row.names=F)





































