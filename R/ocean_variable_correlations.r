
## For calculating correlations among variables shown in Table 2 of:
## "Faulkner, JR, BL Bellerud, DL Widener, SG Smith, and RW Zabel. 2020. Associations among fish length, dam passage history, and survival to adulthood in two at-risk species of Pacific salmon: response to comment."

## This script assumes there is a "data" directory containing the adult return data sets for Chinook and freshwater and ocean data as well as an "output" directory for writing outputs


# the following library needs to be installed and loaded:
library(lme4)


###-----------------------------------
##                 Data
###-----------------------------------

### --- Adult return data -------
cdat.ulgr = read.csv("data/chinook_ULGR_SAR_data.csv", stringsAsFactors=F)
cdat.lgr = read.csv("data/chinook_LGR_SAR_data.csv", stringsAsFactors=F)
# make year a factor
cdat.ulgr$fyear = as.factor(cdat.ulgr$year)
cdat.lgr$fyear = as.factor(cdat.lgr$year)

yrset = 2004:2014
nyears = length(yrset)

## -- Ocean indicators data -- 
# full year set is 1998-2019
odat1 = read.csv("data/ocean_indicators.csv", stringsAsFactors=F)  
# subset to 2004-2014
odat = odat1[odat1$year %in% yrset, ]

### --- Freshwater data -------
# full year set is 1998-2019
idat1 = read.csv("data/annual_freshwater_variables.csv", stringsAsFactors=F)
# subset to 2004-2014
idat = idat1[idat1$year %in% yrset, ]


## --- Merge fresh water and ocean data ---

iodat = merge(odat, idat, all.x=T) 



## --- Summarize returns by year------------
cdat.lgr$ones = 1
cdat.ulgr$ones = 1

# tagged upstream of lgr
mdat.ulgr = aggregate(cdat.ulgr[, c("bypassed", "nbypass", "length", "adult_return") ], by=list(year=cdat.ulgr$fyear), mean)
sdat.ulgr = aggregate(cdat.ulgr[, c("ones", "adult_return") ], by=list(year=cdat.ulgr$fyear), sum)
names(sdat.ulgr)[-1] = c("nfish", "nreturn")
sdat.ulgr$noreturn = sdat.ulgr$nfish - sdat.ulgr$nreturn

# tagged at lgr
mdat.lgr = aggregate(cdat.lgr[, c("bypassed", "nbypass", "length", "adult_return") ], by=list(year=cdat.lgr$fyear), mean)
sdat.lgr = aggregate(cdat.lgr[, c("ones", "adult_return") ], by=list(year=cdat.lgr$fyear), sum)
mdat.lgr$nbypass = mdat.lgr$nbypass + 1
mdat.lgr$bypassed = mdat.lgr$bypassed + 1 
names(sdat.lgr)[-1] = c("nfish", "nreturn")
sdat.lgr$noreturn = sdat.lgr$nfish - sdat.lgr$nreturn

pdat.ulgra = merge(mdat.ulgr, sdat.ulgr, all=T)
pdat.ulgr = merge(pdat.ulgra, iodat, all.x=T)
pdat.ulgr$site = "ulgr"
pdat.lgra = merge(mdat.lgr, sdat.lgr, all=T)
pdat.lgr = merge(pdat.lgra, iodat, all.x=T)
pdat.lgr$site = "lgr"

# combine tagging location data
pdat.comb = rbind(pdat.ulgr, pdat.lgr)
pdat.comb$returns = cbind(pdat.comb$nreturn, pdat.comb$noreturn)

pdat.ulgr$fyear = as.character(pdat.ulgr$year)
pdat.lgr$fyear = as.character(pdat.lgr$year)
pdat.comb$fyear = as.character(pdat.comb$year)

pdat.ulgr$zlength = as.vector(scale(pdat.ulgr$length))
pdat.lgr$zlength = as.vector(scale(pdat.lgr$length))
pdat.comb$zlength = as.vector(scale(pdat.comb$length))

pdat.ulgr$returns = cbind(pdat.ulgr$nreturn, pdat.ulgr$noreturn)
pdat.lgr$returns = cbind(pdat.lgr$nreturn, pdat.lgr$noreturn)

pdat.comb$flow = c(pdat.comb$flow.snk.bon[1:11], pdat.comb$flow.lgr.bon[12:22])
pdat.comb$spillp = c(pdat.comb$spillp.snk.bon[1:11], pdat.comb$spillp.lgr.bon[12:22])
pdat.comb$wtt = c(pdat.comb$wtt.snk.bon[1:11], pdat.comb$wtt.lgr.bon[12:22])

# create standardized data for variables of interest
zdat.comb = scale(pdat.comb[ c(9:24, 27:31, 42:44, 3)])



### -----------------------------------------------
##       Calculate Correlations
### -----------------------------------------------


# calc estimates and pvals for correlations between fresh and salt
cor.est.out = data.frame(variable=dimnames(zdat.comb)[[2]], flow=NA, spillp=NA, wtt=NA, nbypass=NA)
cor.pval.out = data.frame(variable=dimnames(zdat.comb)[[2]], flow=NA, spillp=NA, wtt=NA, nbypass=NA)
for (i in 1:ncol(zdat.comb)) {
  tmp.test = cor.test(zdat.comb[,i], zdat.comb[,"flow"])
  cor.est.out$flow[i] = tmp.test$est
  cor.pval.out$flow[i] = tmp.test$p.val
  tmp.test = cor.test(zdat.comb[,i], zdat.comb[,"spillp"])
  cor.est.out$spillp[i] = tmp.test$est
  cor.pval.out$spillp[i] = tmp.test$p.val
  tmp.test = cor.test(zdat.comb[,i], zdat.comb[,"wtt"])
  cor.est.out$wtt[i] = tmp.test$est
  cor.pval.out$wtt[i] = tmp.test$p.val
  tmp.test = cor.test(zdat.comb[,i], zdat.comb[,"nbypass"])
  cor.est.out$nbypass[i] = tmp.test$est
  cor.pval.out$nbypass[i] = tmp.test$p.val
}

## Write outputs
write.csv(cor.est.out, file="output/ch1_comb_correlation_ocean_fresh.csv", quote=F)
write.csv(cor.pval.out, file="output/ch1_comb_pvals_correlation_ocean_fresh.csv", quote=F)



### -----------------------------------------------
##       Calculate Parameters for Returns
### -----------------------------------------------

## Fit binomial regression models with random year effects and single standardized covariate
## Resulting slope parameters provide estimate of strength of relationship


ysfitout = data.frame(variable=dimnames(zdat.comb)[[2]], estimate=NA, st.err=NA, pval=NA)
for (i in 1:ncol(zdat.comb) ) {
  tmp.fit = glmer(pdat.comb$returns ~ zdat.comb[,i] + (1 | pdat.comb$fyear), family=binomial)
  ysfitout[i, 2:4] = summary(tmp.fit)$coef[2, c(1:2,4)]

}


## write outputs
write.csv(ysfitout, file="output/ch1_comb_survfit_rand_year_ocean_fresh.csv", quote=F)









