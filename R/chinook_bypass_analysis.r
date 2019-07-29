## Chinook Bypass Probability Analysis

## For model fitting in support of:
## Faulkner, JR, BL Bellerud, DL Widener, and RW Zabel. 2019. Associations among fish length, dam passage history, and survival to adulthood in two at-risk species of Pacific salmon.

     
  

library(lme4)



###-------------------------------------------------------------------------------------
##                 Data
###-------------------------------------------------------------------------------------

### --- Read data files
cdat <- read.csv("../../data/chinook_bypass_data.csv", stringsAsFactors=F)


## --- Make new variables
## make year a factor
cdat$fyear <- as.character(cdat$year)




###----------------------------------------------------------------------------------------
##            Model Fitting
###----------------------------------------------------------------------------------------


### ----  Functions ---------------------------------

# function to fit a set of models and output to a list
fitModels <- function(modset, respname, data, subset=NULL){
  nm <- length(modset)
  convec <- numeric(nm)
  outlist <- list()
  for (kk in 1:nm) {
    cat("fitting model: ", kk, " of ", nm, "\n" )
    tmp.form <- formula(paste(respname, "~ ", modset[kk] ))
    tmp.fit <- glmer(tmp.form, family=binomial, data=data, subset=subset, control=glmerControl(optCtrl=list(maxfun=2e5)) )
    if (length(tmp.fit@optinfo$conv$lme4)==0) convec[kk] <- 0
    if (length(tmp.fit@optinfo$conv$lme4)!=0) {
			if (length(grep(x=tmp.fit@optinfo$conv$lme4$messages, pattern="singular"))==1) convec[kk] <- 1
			if (length(grep(x=tmp.fit@optinfo$conv$lme4$messages, pattern="singular"))==0) convec[kk] <- tmp.fit@optinfo$conv$lme4$code
		}
    outlist[[kk]] <- tmp.fit 
    # if did not converge, try to refit from new start
    if (convec[kk]!=0) {
      cat("re-fitting model: ", kk, " of ", nm, "\n" )
      tmp.ss <- getME(tmp.fit,c("theta","fixef"))
      tmp.fit2 <- update(tmp.fit,start=tmp.ss,control=glmerControl(optCtrl=list(maxfun=2e5)) )
      if (length(tmp.fit2@optinfo$conv$lme4)==0) convec[kk] <- 0
	    if (length(tmp.fit2@optinfo$conv$lme4)!=0) {
				if (length(grep(x=tmp.fit2@optinfo$conv$lme4$messages, pattern="singular"))==1) convec[kk] <- 1
				if (length(grep(x=tmp.fit2@optinfo$conv$lme4$messages, pattern="singular"))==0)	convec[kk] <- tmp.fit2@optinfo$conv$lme4$code
			}
      outlist[[kk]] <- tmp.fit2 
    }
    # if still did not converge, try to refit with different optimizer
    if (convec[kk]!=0) {
      cat("re-fitting model: ", kk, " of ", nm, "\n" )
      tmp.ss2 <- getME(tmp.fit2,c("theta","fixef"))
      tmp.fit3 <- update(tmp.fit2,start=tmp.ss2,control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
      if (length(tmp.fit3@optinfo$conv$lme4)==0) convec[kk] <- 0
	    if (length(tmp.fit3@optinfo$conv$lme4)!=0) {
				if (length(grep(x=tmp.fit3@optinfo$conv$lme4$messages, pattern="singular"))==1) convec[kk] <- 1
				if (length(grep(x=tmp.fit3@optinfo$conv$lme4$messages, pattern="singular"))==0)	convec[kk] <- tmp.fit3@optinfo$conv$lme4$code
			}
      outlist[[kk]] <- tmp.fit3 
    }
    if (convec[kk]!=0) cat("2nd re-fit of model ", kk, " failed!", "\n" )
  }
 outlist$converge <- convec
 return(outlist)
}



## Function to create table of AIC ranks and other stats
getAICsummary <- function(fit.list, mod.list, subset=NULL) {
	nmod <- length(fit.list)-1
	modnums <- 1:nmod
	if (!is.null(subset)) {
		modnums <- modnums[subset]
		nmod <- length(modnums)
		mod.list <- mod.list[subset]
	}
	aicv1 <- as.vector((sapply(fit.list[modnums], FUN=extractAIC))[2,])   #extract aic
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
# S = release site 
# G = day at Lower Granite 
# G2 = G squared
# D = day at Bonneville or Trawl
# D2 = D squared
# L = length


### Possible models models
# notes: - all models have Y and C and S 
#        - G variables only occur in Snake River models and D variables only occur in Columbia River models
#        - All models are constructed separately by Dam


# Snake models
# 1) Y + C + S
# 2) Y + C + S + G
# 3) Y + C + S + G + G2
# 4) Y + C + S + G + Y*G
# 5) Y + C + S + G + G2 + Y*G
# 6) Y + C + S + G + G2 + Y*G + Y*G2



# Columbia models
# 2) Y + C + S
# 5) Y + C + S + D
# 6) Y + C + S + D + D2
# 10) Y + C + S + D + Y*D
# 11) Y + C + S + D + D2 + Y*D
# 12) Y + C + S + D + D2 + Y*D + Y*D2




### --- Models Without Length

## Snake models
modsets.sr1 <- c("clipped + rel_site2 + (1 | fyear)",
"clipped + rel_site2 + zday_lgr + (1 | fyear)",
"clipped + rel_site2 + zday_lgr + zday_lgr_sq + (1 | fyear)", 
"clipped + rel_site2 + zday_lgr + (zday_lgr | fyear)",
"clipped + rel_site2 + zday_lgr + zday_lgr_sq + (zday_lgr | fyear)",
"clipped + rel_site2 + zday_lgr + zday_lgr_sq + (zday_lgr + zday_lgr_sq | fyear)")


# Columbia models
modsets.cr1 <- c("clipped + rel_site2 + (1 | fyear)",
"clipped + rel_site2 + zlast_day + (1 | fyear)",
"clipped + rel_site2 + zlast_day + zlast_day_sq + (1 | fyear)", 
"clipped + rel_site2 + zlast_day + (zlast_day | fyear)",
"clipped + rel_site2 + zlast_day + zlast_day_sq + (zlast_day | fyear)",
"clipped + rel_site2 + zlast_day + zlast_day_sq + (zlast_day + zlast_day_sq | fyear)" )


### --- Models With Length
modsets.sr <- paste("zlength +", modsets.sr1)
modsets.cr <- paste("zlength +", modsets.cr1)



nmods <- length(modsets.sr)
modnums <- 1:nmods

respsets <- c("bypass_LGR", "bypass_LGS", "bypass_LMN", "bypass_IHR", "bypass_MCN", "bypass_JDA", "bypass_BON", "bypass_BON")

damsets <- c("LGR", "LGS", "LMN", "IHR", "MCN", "JDA", "BON", "BCC")

subsets <- list(cdat$rel_site2!="LGR", cdat$last_site=="twx", (cdat$bypass_BON==1 | cdat$pass_BCC==1) & cdat$fyear %in% as.character(2006:2014))





###--------------------------------------------------
##        Lower Granite Dam 
###--------------------------------------------------

idam <- which(damsets=="LGR")

# -------   Fit models with length -------------------------------------
modfits.lgr <- fitModels(modsets.sr, respsets[idam], cdat, subsets[[1]])

# check if any models did not converge (any non-zero)
modfits.lgr$converge 
# all are fine


# Calculate AIC
aicsum.lgr <- getAICsummary(modfits.lgr, modsets.sr, subset = modfits.lgr$converge==0 )
aicsum.lgr


# find best model
bestmod.lgr <- aicsum.lgr$Model[1]
modfit.best.lgr <- modfits.lgr[[bestmod.lgr]]
modsum.best.lgr <- summary(modfit.best.lgr)
modsum.best.lgr


# -------   Fit best model without length -------------------------------------
tmpset <- modsets.sr1[bestmod.lgr]
modfit.nolen.lgr <- fitModels(tmpset, respsets[idam], cdat, subsets[[1]])
modfit.nolen.lgr$converge
modsum.nolen.lgr <- summary(modfit.nolen.lgr[[1]])
modsum.nolen.lgr




###--------------------------------------------------
##        Little Goose Dam 
###--------------------------------------------------

idam <- which(damsets=="LGS")

# -------   Fit models with length -------------------------------------
modfits.lgs <- fitModels(modsets.sr, respsets[idam], cdat)

# check if any models did not converge (any non-zero)
modfits.lgs$converge 
# all are fine

# Calculate AIC
aicsum.lgs <- getAICsummary(modfits.lgs, modsets.sr, subset = modfits.lgs$converge==0 )
aicsum.lgs


# find best model
bestmod.lgs <- aicsum.lgs$Model[1]
modfit.best.lgs <- modfits.lgs[[bestmod.lgs]]
modsum.best.lgs <- summary(modfit.best.lgs)
modsum.best.lgs

# -------   Fit best model without length -------------------------------------
tmpset <- modsets.sr1[bestmod.lgs]
modfit.nolen.lgs <- fitModels(tmpset, respsets[idam], cdat)
modfit.nolen.lgs$converge
modsum.nolen.lgs <- summary(modfit.nolen.lgs[[1]])
modsum.nolen.lgs






###--------------------------------------------------
##        Lower Monumental Dam 
###--------------------------------------------------

idam <- which(damsets=="LMN")

# -------   Fit models with length -------------------------------------
modfits.lmn <- fitModels(modsets.sr, respsets[idam], cdat)

# check if any models did not converge (any non-zero)
modfits.lmn$converge
# all are fine

# Calculate AIC
aicsum.lmn <- getAICsummary(modfits.lmn, modsets.sr, subset = modfits.lmn$converge==0 )
aicsum.lmn


# find best model
bestmod.lmn <- aicsum.lmn$Model[1]
modfit.best.lmn <- modfits.lmn[[bestmod.lmn]]
modsum.best.lmn <- summary(modfit.best.lmn)
modsum.best.lmn

# -------   Fit best model without length -------------------------------------
tmpset <- modsets.sr1[bestmod.lmn]
modfit.nolen.lmn <- fitModels(tmpset, respsets[idam], cdat)
modfit.nolen.lmn$converge
modsum.nolen.lmn <- summary(modfit.nolen.lmn[[1]])
modsum.nolen.lmn




###--------------------------------------------------
##        Ice Harbor Dam 
###--------------------------------------------------

idam <- which(damsets=="IHR")

# -------   Fit models with length -------------------------------------
modfits.ihr <- fitModels(modsets.sr, respsets[idam], cdat, cdat$fyear %in% as.character(2005:2014))

# check if any models did not converge (any non-zero)
modfits.ihr$converge 

# Calculate AIC
aicsum.ihr <- getAICsummary(modfits.ihr, modsets.sr, subset = modfits.ihr$converge==0 )
aicsum.ihr


# find best model
bestmod.ihr <- aicsum.ihr$Model[1]
modfit.best.ihr <- modfits.ihr[[bestmod.ihr]]
modsum.best.ihr <- summary(modfit.best.ihr)
modsum.best.ihr


# -------   Fit best model without length -------------------------------------
tmpset <- modsets.sr1[bestmod.ihr]
modfit.nolen.ihr <- fitModels(tmpset, respsets[idam], cdat, cdat$fyear %in% as.character(2005:2014))
modfit.nolen.ihr$converge
modsum.nolen.ihr <- summary(modfit.nolen.ihr[[1]])
modsum.nolen.ihr




###--------------------------------------------------
##        McNary Dam 
###--------------------------------------------------

idam <- which(damsets=="MCN")

# -------   Fit models with length -------------------------------------
modfits.mcn <- fitModels(modsets.cr, respsets[idam], cdat)

# check if any models did not converge (any non-zero)
modfits.mcn$converge 
# all are fine

# Calculate AIC
aicsum.mcn <- getAICsummary(modfits.mcn, modsets.cr, subset = modfits.mcn$converge==0 )
aicsum.mcn


# find best model
bestmod.mcn <- aicsum.mcn$Model[1]
modfit.best.mcn <- modfits.mcn[[bestmod.mcn]]
modsum.best.mcn <- summary(modfit.best.mcn)
modsum.best.mcn


# -------   Fit best model without length -------------------------------------
tmpset <- modsets.cr1[bestmod.mcn]
modfit.nolen.mcn <- fitModels(tmpset, respsets[idam], cdat)
modfit.nolen.mcn$converge

modsum.nolen.mcn <- summary(modfit.nolen.mcn[[1]])
modsum.nolen.mcn




###--------------------------------------------------
##        John Day Dam 
###--------------------------------------------------

idam <- which(damsets=="JDA")

# -------   Fit models with length -------------------------------------
modfits.jda <- fitModels(modsets.cr, respsets[idam], cdat)

# check if any models did not converge (any non-zero)
modfits.jda$converge 
# all are fine

# Calculate AIC
aicsum.jda <- getAICsummary(modfits.jda, modsets.cr, subset = modfits.jda$converge==0 )
aicsum.jda


# find best model
bestmod.jda <- aicsum.jda$Model[1]
modfit.best.jda <- modfits.jda[[bestmod.jda]]
modsum.best.jda <- summary(modfit.best.jda)
modsum.best.jda


# -------   Fit best model without length -------------------------------------
tmpset <- modsets.cr1[bestmod.jda]
modfit.nolen.jda <- fitModels(tmpset, respsets[idam], cdat)
modfit.nolen.jda$converge

modsum.nolen.jda <- summary(modfit.nolen.jda[[1]])
modsum.nolen.jda



###--------------------------------------------------
##        Bonneville Dam 
###--------------------------------------------------

idam <- which(damsets=="BON")

# -------   Fit models with length -------------------------------------
modfits.bon <- fitModels(modsets.cr, respsets[idam], cdat, subsets[[2]] )

# check if any models did not converge (any non-zero)
modfits.bon$converge 

# Calculate AIC
aicsum.bon <- getAICsummary(modfits.bon, modsets.cr, subset = modfits.bon$converge==0 )
aicsum.bon


# find best model
bestmod.bon <- aicsum.bon$Model[1]
modfit.best.bon <- modfits.bon[[bestmod.bon]]
modsum.best.bon <- summary(modfit.best.bon)
modsum.best.bon


# -------   Fit best model without length -------------------------------------
tmpset <- modsets.cr1[bestmod.bon]
modfit.nolen.bon <- fitModels(tmpset, respsets[idam], cdat, subsets[[2]])
modfit.nolen.bon$converge
modsum.nolen.bon <- summary(modfit.nolen.bon[[1]])
modsum.nolen.bon




###--------------------------------------------------
##        Bonneville Dam Corner Collector
###--------------------------------------------------


idam <- which(damsets=="BCC")

# -------   Fit models with length -------------------------------------
modfits.bcc <- fitModels(modsets.cr, respsets[idam], cdat, subsets[[3]])

# check if any models did not converge (any non-zero)
modfits.bcc$converge 
# all are fine

# Calculate AIC
aicsum.bcc <- getAICsummary(modfits.bcc, modsets.cr, subset = modfits.bcc$converge==0 )
aicsum.bcc


# find best model 
bestmod.bcc <- aicsum.bcc$Model[1]
modfit.best.bcc <- modfits.bcc[[bestmod.bcc]]
modsum.best.bcc <- summary(modfit.best.bcc)
modsum.best.bcc

# -------   Fit best model without length -------------------------------------
tmpset <- modsets.cr1[bestmod.bcc]
modfit.nolen.bcc <- fitModels(tmpset, respsets[idam], cdat, subsets[[3]] )
modfit.nolen.bcc$converge
modsum.nolen.bcc <- summary(modfit.nolen.bcc[[1]])
modsum.nolen.bcc




## -----------------------------------------
##     Create Summary for Table 1
## -----------------------------------------

# species, reartype, dam, n, estimate, 95%CI, dAIC, p-value

dnames <- c("LGD","LGSD","LMD","IHD","MCD","JDD","BVD","BVCC")

tab1 <- data.frame(Species=rep("Chinook",8), Dam=dnames, 
  n = c(length(modsum.best.lgr$resid),
    length(modsum.best.lgs$resid),
    length(modsum.best.lmn$resid),
    length(modsum.best.ihr$resid),
    length(modsum.best.mcn$resid),
    length(modsum.best.jda$resid),
    length(modsum.best.bon$resid),
    length(modsum.best.bcc$resid) ),
  round(rbind(extractCI(modsum.best.lgr, "zlength", scale="log.odds"),
    extractCI(modsum.best.lgs, "zlength", scale="log.odds"),
    extractCI(modsum.best.lmn, "zlength", scale="log.odds"),
    extractCI(modsum.best.ihr, "zlength", scale="log.odds"),
    extractCI(modsum.best.mcn, "zlength", scale="log.odds"),
    extractCI(modsum.best.jda, "zlength", scale="log.odds"),
    extractCI(modsum.best.bon, "zlength", scale="log.odds"),
    extractCI(modsum.best.bcc, "zlength", scale="log.odds") ),3 ),
  dAIC=c(round(modsum.best.lgr$AIC[1] - modsum.nolen.lgr$AIC[1], 1),
    round(modsum.best.lgs$AIC[1] - modsum.nolen.lgs$AIC[1], 1),
    round(modsum.best.lmn$AIC[1] - modsum.nolen.lmn$AIC[1], 1),
    round(modsum.best.ihr$AIC[1] - modsum.nolen.ihr$AIC[1], 1),
    round(modsum.best.mcn$AIC[1] - modsum.nolen.mcn$AIC[1], 1),
    round(modsum.best.jda$AIC[1] - modsum.nolen.jda$AIC[1], 1),
    round(modsum.best.bon$AIC[1] - modsum.nolen.bon$AIC[1], 1),
    round(modsum.best.bcc$AIC[1] - modsum.nolen.bcc$AIC[1], 1) ) ,
  pvalue=c(round(modsum.best.lgr$coef["zlength", 4], 4),
    round(modsum.best.lgs$coef["zlength", 4], 4),
    round(modsum.best.lmn$coef["zlength", 4], 4),
    round(modsum.best.ihr$coef["zlength", 4], 4),
    round(modsum.best.mcn$coef["zlength", 4], 4),
    round(modsum.best.jda$coef["zlength", 4], 4),
    round(modsum.best.bon$coef["zlength", 4], 4),
    round(modsum.best.bcc$coef["zlength", 4], 4) ) 
)

names(tab1)[4:6] <- c("Estimate", "95%CI.L", "95%CI.U")

tab1


write.csv(tab1, file="output/chinook_Table1_values.csv", quote=F, row.names=F)





## ------------------------------------------
##     Create Summary for Table A3 
## ------------------------------------------

tabA3 <- data.frame(Species=rep("Chinook",8),  Dam=dnames, 
  modform=c(aicsum.lgr$ModelForm[1],
    aicsum.lgs$ModelForm[1],
    aicsum.lmn$ModelForm[1],
    aicsum.ihr$ModelForm[1],
    aicsum.mcn$ModelForm[1],
    aicsum.jda$ModelForm[1],
    aicsum.bon$ModelForm[1],
    aicsum.bcc$ModelForm[1])
)

tabA3

write.csv(tabA3, file="output/chinook_TableA3_values.csv", quote=F, row.names=F)






