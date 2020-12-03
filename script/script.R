# load libraries
library(tidyverse)
library(mrgsolve)

# Chunks 1-2: simplePBPK
# Chunks 3-10: voriPBPK model
# Chunks 11-12: sensitivity analysis
# Chunks 13-14: parameter estimation
# Chunks 15-18: population simulation
# Chunk 19: interactive presentation / shiny


#------------------------------------------------------------------------------#

### Simple PBPK ###

################################################################################
################################## Chunk 1  ####################################
################################################################################

# Populate and compile simplePBPK model
mod <- mread("simplePBPK", "../model")

################################################################################
################################################################################

################################################################################
################################## Chunk 2 #####################################
################################################################################

## Simulate a 100 mg dose given as an IV bolus
mod %>%
  ev(amt=100, cmt="VEN") %>%
  mrgsim(end=100) %>%
  plot()

################################################################################
################################################################################


#------------------------------------------------------------------------------#

### voriconazole PBPK ###

################################################################################
################################## Chunk 3 #####################################
################################################################################

# Compile voriPBPK model
modA <- mread("../model/voriPBPK")

################################################################################
################################################################################

################################################################################
##################################  Chunk 4 ####################################
################################################################################

# Use `calcKp_PT.R` function to calculate voriconzole tissue:plasma partition 
# coefficients according to Poulin and Theil method https://www.ncbi.nlm.nih.gov/pubmed/11782904.

# source partition coefficient calculation script
source("calcKp_PT.R")

#voriconazole physicochemical properties
logP <- 2.56  #lipophilicity
pKa <- 1.76  
fup <- 0.42   #unbound fraction in plasma
type <- 3     #monoprotic base
BP <- 1       #blood:plasma concentration ratio
dat <- read.csv("../data/source/tissue_comp_P&T.csv")

#calculate partition coefficients
Kp <- calcKp_PT(logP=logP, pKa=pKa, fup=fup, BP=BP, type=type, dat=dat)

#update model parameters partition coefficients
modA <- param(modA, Kp)

################################################################################
################################################################################

################################################################################
##################################  Chunk 5 ####################################
################################################################################

# Simulate the steady state after a 4 mg/kg IV infusion dose of voriconazole given 
# to an adult male with a rate of 4 mg/kg/h twice a day for 7 days. Compare the 
# steady state plasma drug concentration-time profiles from previous simulation to 
# the observed data in `Adult_IV.csv`. (N.B.: observed data were digitized from 
# Zane and Thakker (2014) paper using WebPlotDigitizer https://automeris.io/WebPlotDigitizer/):
  
#load observed IV infusion data
obs <- read.csv("../data/source/Adult_IV.csv")

#set simulation conditions
bw   <- 73
amt  <- 4*bw
rate <- 4*bw
cmt  <- "VEN"
ii   <- 12
addl <- 13
ss   <- 1

#run simulation
sim <- 
  modA %>% 
  ev(amt=amt, cmt=cmt, ii=ii, addl=addl, rate=rate, ss=ss) %>% 
  mrgsim(delta = 0.1, end = 12) %>% 
  filter(row_number() != 1)  

#plot prediction and compare to observed data
gp <- ggplot() + 
  geom_point(data = obs, aes(x=time, y=obs, col="observed"), size=2.5) + 
  geom_errorbar(data = obs, aes(x = time, y = obs, ymin=obs-sd, ymax=obs+sd), width=0) +
  geom_line(data = sim, aes(x=time, y=CP, col="sim"), lwd=1) + 
  scale_colour_manual(name='', 
                      values=c('sim'='black', 'observed'='black'), 
                      breaks=c("observed","sim"),
                      labels=c("observed","predicted")) +
  guides(colour = guide_legend(override.aes = list(linetype=c(0,1), shape=c(16, NA)))) +
  labs(title="Adult 4 mg/kg IV", x="time (h)", y="Plasma concentration (mg/L)") +
  theme_bw()
gp

################################################################################
################################################################################

################################################################################
##################################  Chunk 6 ####################################
################################################################################

# simulate the steady state after a 200 mg PO dose of voriconazole given to an adult male 
# with a rate of twice a day for 7 days. Compare the steady state plasma drug 
# concentration-time profile to the observed data in `Adult_PO.csv`.

obs <- read.csv("../data/source/Adult_PO.csv")

bw   <- 73
amt  <- 200
cmt  <- "GUTLUMEN"
ii   <- 12
addl <- 13
ss   <- 1

sim <- 
  modA %>% 
  ev(amt=amt, cmt=cmt, ii=ii, addl=addl, ss=ss) %>% 
  mrgsim(delta = 0.1, end = 12) %>% 
  filter(row_number() != 1)  

gp <- ggplot() + 
  geom_point(data = obs, aes(x=time, y=obs, col="observed"), size=2.5) + 
  geom_errorbar(data = obs, aes(x = time, y = obs, ymin=obs-sd, ymax=obs+sd), width=0) +
  geom_line(data = sim, aes(x=time, y=CP, col="sim"), lwd=1) + 
  scale_colour_manual(name='', 
                      values=c('sim'='black', 'observed'='black'), 
                      breaks=c("observed","sim"),
                      labels=c("observed","predicted")) +
  guides(colour = guide_legend(override.aes = list(linetype=c(0,1), shape=c(16, NA)))) +
  labs(title="Adult 4 mg/kg PO", x="time (h)", y="Plasma concentration (mg/L)") +
  theme_bw()
gp

################################################################################
################################################################################

################################################################################
##################################  Chunk 7 ####################################
################################################################################

# Integrate the model with an additional gut wall enterocyte compartment to account for 
# intestinal clearance, intestinal transit and lumen solubility effects on absorption rate. 
# Note: use about an intestinal clearance that is 30 times lower than hepatic clearance.
# Recompile and re-run the previous step. Any change!!
# https://www.jstage.jst.go.jp/article/yakushi/123/5/123_5_369/_article/-char/en

modA <- mread("../model/voriPBPK_ext") %>%
  param(MPPGI = 30.3/30) %>%
  param(Kp)

sim <- 
  modA %>% 
  ev(amt=amt, cmt=cmt, ii=ii, addl=addl, ss=ss) %>% 
  mrgsim(delta = 0.1, end = 12) %>% 
  filter(row_number() != 1)  

gp <- ggplot() + 
  geom_point(data = obs, aes(x=time, y=obs, col="observed"), size=2.5) + 
  geom_errorbar(data = obs, aes(x = time, y = obs, ymin=obs-sd, ymax=obs+sd), width=0) +
  geom_line(data = sim, aes(x=time, y=CP, col="sim"), lwd=1) + 
  scale_colour_manual(name='', 
                      values=c('sim'='black', 'observed'='black'), 
                      breaks=c("observed","sim"),
                      labels=c("observed","predicted")) +
  guides(colour = guide_legend(override.aes = list(linetype=c(0,1), shape=c(16, NA)))) +
  labs(title="Adult 4 mg/kg PO", x="time (h)", y="Plasma concentration (mg/L)") +
  theme_bw()
gp

################################################################################
################################################################################


################################################################################
##################################  Chunk 8 ####################################
################################################################################

# Generate pediatric model

# pediatric (5 yo) male physiology; https://www.ncbi.nlm.nih.gov/pubmed/14506981
pedPhys <- list(WEIGHT = 19,
                Vad = 5.5,
                Vbo = 2.43,
                Vbr = 1.31,
                VguWall = 0.22,
                VguLumen = 0.117,
                Vhe = 0.085,
                Vki = 0.11,
                Vli = 0.467,
                Vlu = 0.125,
                Vmu = 5.6,
                Vsp = 0.05,
                Vbl = 1.5,
                Qad = 0.05*3.4*60,
                Qbo = 0.05*3.4*60,
                Qbr = 0.12*3.4*60,
                Qgu = 0.15*3.4*60, 
                Qhe = 0.04*3.4*60,
                Qki = 0.19*3.4*60,
                Qmu = 0.17*3.4*60,
                Qsp = 0.03*3.4*60,
                Qha = 0.065*3.4*60, 
                Qlu = 3.4*60,
                MPPGL = 26,
                VmaxH = 120.5,
                KmH = 11,
                MPPGI = 0,
                VmaxG = 120.5,
                KmG = 11)

modP <- param(modA, pedPhys)

################################################################################
################################################################################

################################################################################
##################################  Chunk 9 ####################################
################################################################################

# Simulate the steady state after a 4 mg/kg voriconazole IV infusion dosing in a male child 
# subject infused at a rate of 3 mg/kg/h twice a day for seven days. Compare the steady state 
# to the observed data in `Pediatric_IV.csv`.

obs <- read.csv("../data/source/Pediatric_IV.csv")  #load observed data

wt   <- 19  #pediatric body weight
amt  <- 4*wt  
rate <- 3*wt
cmt  <- "VEN"  #intravenous infusion
ii   <- 12
addl <- 13
ss   <- 1

# simulate
sim <- 
  modP %>%
  ev(amt=amt, cmt=cmt, ii=ii, addl=addl, rate=rate, ss=1) %>% 
  mrgsim(delta = 0.1, end = 12) %>% 
  dplyr::filter(row_number() != 1)  

# plot
gp <- ggplot() + 
  geom_point(data = obs, aes(x=time, y=obs, col="observed"), size=2.5) + 
  geom_errorbar(data = obs, aes(x = time, y = obs, ymin=obs-sd, ymax=obs+sd), width=0) +
  geom_line(data = sim, aes(x=time, y=CP, col="sim"), lwd=1) + 
  scale_colour_manual(name='', 
                      values=c('sim'='black', 'observed'='black'), 
                      breaks=c("observed","sim"),
                      labels=c("observed","predicted")) +
  guides(colour = guide_legend(override.aes = list(linetype=c(0,1), shape=c(16, NA)))) +
  labs(title="Pediatric 4 mg/kg IV", x="time (h)", y="Plasma concentration (mg/L)") +
  theme_bw()
gp

################################################################################
################################################################################

################################################################################
#################################  Chunk 10 ####################################
################################################################################

# Simulate the steady state after a 4 mg/kg voriconazole PO dosing in a male child subject 
# twice a day for seven days. Compare to obsreved data in `Pediatric_PO.csv`
# Note: Include a similar 30-fold lower intestinal clearance than hepatic clearance.
# https://www.jstage.jst.go.jp/article/yakushi/123/5/123_5_369/_article/-char/en

obs <- read.csv("../data/source/Pediatric_PO.csv")  #load observed data

# adjust intestinal clearance
modP <- modP %>% param(MPPGI = 26 / 30)

# simulation conditions
bw   <- 19
amt  <- 4 * bw
cmt  <- "GUTLUMEN"
ii   <- 12
addl <- 13
ss   <- 1

# simulate
sim <- 
  modP %>%
  ev(amt=amt, cmt=cmt, ii=ii, addl=addl, ss=1) %>% 
  mrgsim(delta = 0.1, end = 12) %>% 
  dplyr::filter(row_number() != 1)  

# plot
gp <- ggplot() + 
  geom_point(data = obs, aes(x=time, y=obs, col="observed"), size=2.5) + 
  geom_errorbar(data = obs, aes(x = time, y = obs, ymin=obs-sd, ymax=obs+sd), width=0) +
  geom_line(data = sim, aes(x=time, y=CP, col="sim"), lwd=1) + 
  scale_colour_manual(name='', 
                      values=c('sim'='black', 'observed'='black'), 
                      breaks=c("observed","sim"),
                      labels=c("observed","predicted")) +
  guides(colour = guide_legend(override.aes = list(linetype=c(0,1), shape=c(16, NA)))) +
  labs(title="Pediatric 4 mg/kg PO", x="time (h)", y="Plasma concentration (mg/L)") +
  theme_bw()
gp

################################################################################
################################################################################

#------------------------------------------------------------------------------#

### Sensitivity Analysis ###

################################################################################
#################################  Chunk 11 ####################################
################################################################################

## graphical ##

# Run graphical sensitivity analysis for the muscle:plasma (`Kpmu`) and 
# lung:plasma (`Kplu`) partition coefficients using adult IV data.

#load observed IV infusion data
obs <- read.csv("../data/source/Adult_IV.csv")

#set simulation conditions
bw   <- 73
amt  <- 4*bw
rate <- 4*bw
cmt  <- "VEN"
ii   <- 12
addl <- 13
ss   <- 1

##' Define an intervention
e <- ev(cmt=cmt, amt=amt, rate=rate, ii= ii, addl=addl, ss=1)

## Sensitivity analysis on Kpmu
idata <- expand.idata(Kpmu = c(3/2, 3, 3*2))

modA %>%
  carry_out(Kpmu) %>% 
  mrgsim_ei(e, idata, delta = 0.1, recsort=3, obsonly=TRUE, end = 12) %>% 
  mutate(Kpmu = factor(Kpmu)) %>%
  ggplot(aes(x=time, y=CP, col=Kpmu)) +
  geom_line() +
  labs(x="time (h)", y="Plasma concentration (mg/L)") +
  theme_bw()

##' Sensitivity analysis on Kplu
idata <- expand.idata(Kplu = c(1/2, 1, 1*2))

modA %>%
  carry_out(Kplu) %>% 
  mrgsim_ei(e, idata, delta = 0.1, recsort=3, obsonly=TRUE, end = 12) %>% 
  mutate(Kplu = factor(Kplu)) %>%
  ggplot(aes(x=time, y=CP, col=Kplu)) +
  geom_line() +
  labs(x="time (h)", y="Plasma concentration (mg/L)") +
  theme_bw()

################################################################################
################################################################################

################################################################################
#################################  Chunk 12 ####################################
################################################################################

## FME package ##

# Use package `FME` (https://cran.r-project.org/web/packages/FME/index.html) to 
# run local sensitivity analysis for the partition coefficient parameters Soetaert, 
# Karline, and Thomas Petzoldt. 2010. “Inverse Modelling, Sensitivity and Monte Carlo 
# Analysis in R Using Package FME.” Journal of Statistical Software, Articles 33 (3): 1–28. 
# https://cran.r-project.org/web/packages/FME/vignettes/FME.pdf. <br> What is the 
# most influential parameter?

library(FME)

### sensitivity analysis
#set a list of the candidate parameters
parms <- list(Kpad = 9.89, Kpbo = 7.91, Kpbr = 7.35, Kpgu = 5.82, Kphe = 1.95, Kpki = 2.9, Kpli = 4.66, Kplu = 0.83, Kpmu = 2.94, Kpsp = 2.96) 

#set the output variable of interest
sensvar <- c("CP")  

#set the output function
outFun <- function(pars){
  out <- modA %>%
    param(pars) %>%
    ev(amt=amt, cmt=cmt, ii=ii, addl=addl, rate=rate, ss=ss) %>%
    mrgsim(end=12, delta=0.1) %>%
    dplyr::filter(row_number() != 1) 
  
  #get rid of ID because sensFun will take the first column as the time variable
  out <- out %>% dplyr::select(-ID)  
  return(out)
}

locSens <- sensFun(func=outFun, parms=parms, sensvar=sensvar, tiny=1e-4)
summary(locSens)

## summary plots
plot(summary(locSens))

# nicer view
summ <- as_tibble(summary(locSens)) %>%
  mutate(parms = names(parms)) 

ggplot(data=summ, aes(x=reorder(parms, Mean), y=Mean)) + 
  geom_col() + 
  labs(x="Parameter", y="Coefficient") +
  coord_flip() +
  geom_hline(yintercept = 0, lty=2) +
  theme_bw()

## time-course sensitivity
plot(locSens, which = c("CP"), xlab="time", lwd = 2)

#nicer view
df_temp <- as_tibble(locSens) %>%
  gather(Parameter, Coefficient, -x, -var) %>%
  mutate(Parameter = factor(Parameter)) %>%
  rename(time=x) %>%
  group_by(Parameter) %>%
  mutate(Coefficient = Coefficient - first(Coefficient)) %>%
  ungroup()

ggplot(data=df_temp, aes(x=time, y=Coefficient, col=Parameter)) +
  geom_line() +
  theme(legend.position="right") +
  facet_wrap(~var) +
  theme_bw()

################################################################################
################################################################################

#------------------------------------------------------------------------------#

### Parameter estimation ###

################################################################################
################################## Chunk 13 ####################################
################################################################################

# Use package `nloptr` (https://cran.r-project.org/web/packages/nloptr/index.html) 
# to optimize for most influential partition coefficient parameter and compare prediction 
# before and after optimization to observed data
# Use package `numDeriv` (https://cran.r-project.org/web/packages/numDeriv/index.html) 
# to generate the 95% CI around the parameter estimates
# For more background on optimization methods, check out Metrum's open science course 
# MI210 https://metrumrg.com/course/mi210-essentials-population-pk-pd-modeling-simulation/

library(nloptr)
library(numDeriv)
library(kableExtra)

### optimization with maximum likelihood method
theta <- log(c(Kpmu=2.94, sigma=1))  #initial parameter for MLE; mean and variance
sampl <- obs$time                    #sampling times from observed data

bw <- 73
amt = 4*bw
rate = 4*bw
cmt = "VEN"
ii = 12
addl = 13
ss = 1

##set up objective function
OF <- function(pars, pred=F){
  pars <- lapply(pars,exp)  #Get out of log domain for MLE
  pars <- as.list(pars)
  names(pars) <- names(theta)
  
  ## Get a prediction
  out <- as.data.frame(modA %>% 
                         param(pars) %>%
                         ev(cmt=cmt, amt=amt, rate=rate, ii=ii, addl=addl, ss=ss) %>%
                         mrgsim(end=-1, add=sampl)
  )
  
  out <- out[-1,] 
  
  if(pred) return(out)
  
  ##OLS
  #return(sum((out$Cvenous - obs$obs)^2))
  
  ##maximum likelihood
  return(-1*sum(dnorm(log(obs$obs),
                      mean=log(out$CP),
                      sd=pars$sigma, log=TRUE)))
}

##Fit with nloptr package
##derivative-free optimizers
#fit <- neldermead(theta, OF)  #Nelder-Mead simplex
fit <- newuoa(theta, OF)  #New Unconstrained Optimization with quadratic Approximation

##gradient-basd optimizers
#fit <- tnewton(theta, OF)  #Local optimizer; Nelder-Mead simplex
#fit <- mlsl(theta, OF, lower=log(c(0.1)), upper=log(c(10)))  #global optimizer; multi-level single-linkage; takes very long ~ 10 min

##global optimizers
#fit <- direct(OF, lower=log(c(0.1, 0.1)), upper=log(c(10, 2)))  #DIviding RECTangles algorithm; takes ~ 3 min

fit

p <- as.list(exp(fit$par))  #get the parameters on the linear scale
names(p) <- names(theta)
p

# get standard error and confidence intervals around the estimated parameters
h <- hessian(OF, fit$par)
vc_log <- solve(h)  #variance-covariance matrix
SE_log <- sqrt(diag(vc_log))  #standard error on log scale

# create dataframe with parameter summary
sig <- function(x) signif(x, 3)
paramSumm <- tibble(Parameter = names(theta),
                    Estimate = exp(fit$par),
                    lb = exp(fit$par - (1.96 * SE_log)),
                    ub = exp(fit$par + (1.96 * SE_log))) %>%
  mutate(`Estimate (95% CI)` = paste0(sig(Estimate), " (", sig(lb), ", ", sig(ub), ")")) %>%
  select(Parameter, `Estimate (95% CI)`)

paramSumm %>%
  knitr::kable() %>%
  kable_styling()

# compare initial predictions to those using the optimized parameters
predB4 <- OF(theta, pred=T)
predAfter <- OF(fit$par, pred=T)

gp <- ggplot() +
  geom_point(data=obs, aes(x=time, y=obs)) + 
  geom_errorbar(data=obs, aes(x=time, y=obs, ymin=obs-sd, ymax=obs+sd), width=0) +
  geom_line(data=predB4, aes(x=time, y=CP), lty=2) +
  geom_line(data=predAfter, aes(x=time, y=CP)) +
  labs(x="time (h)", y="Plasma concentration (mg/L)") + 
  theme_bw()
gp


################################################################################
################################################################################

################################################################################
################################## Chunk 14 ####################################
################################################################################

# Simulate a 4 mg/kg voriconazole IV infusion dosing in a male child subject infused 
# at a rate of 3 mg/kg/h twice a day for seven days. Compare the steady state prediction 
# before and after optimization to the observed data in `Pediatric_IV.csv`.

obs <- read.csv("../data/source/Pediatric_IV.csv")  #load observed data

wt <- 19  #adult body weight
amt <- 4*wt  
rate <- 3*wt
cmt <- "VEN"  #intravenous infusion
ii = 12
addl = 13
ss = 1

# simulate
simB4 <- as.data.frame(modP %>%
                         ev(amt=amt, cmt=cmt, ii=ii, addl=addl, rate=rate, ss=1) %>% 
                         mrgsim(delta = 0.1, end = 12)) %>% 
  dplyr::filter(row_number() != 1)  

simAfter <- as.data.frame(modP %>%
                            param(Kpmu=p$Kpmu) %>%
                            ev(amt=amt, cmt=cmt, ii=ii, addl=addl, rate=rate, ss=1) %>% 
                            mrgsim(delta = 0.1, end = 12)) %>% 
  dplyr::filter(row_number() != 1)

# plot
gp <- ggplot() + 
  geom_point(data = obs, aes(x=time, y=obs), size=2.5) + 
  geom_errorbar(data = obs, aes(x = time, y = obs, ymin=obs-sd, ymax=obs+sd), width=0) +
  geom_line(data = simB4, aes(x=time, y=CP), lty=2) + 
  geom_line(data = simAfter, aes(x=time, y=CP)) + 
  labs(title="Pediatric 4 mg/kg IV", x="time (h)", y="Plasma concentration (mg/L)") +
  theme_bw()
gp

################################################################################
################################################################################

#------------------------------------------------------------------------------#

### Population simulation ###

################################################################################
################################## Chunk 15 ####################################
################################################################################

# Simulate a simple voriconazole dosing scenario for a population of 100 individuals, 
# 50 males and 50 females, with ages ranging between 20-80, weights ranging between 
# 50-100 kg and heights between 1.5 and 1.9 m. Use population saved in file 
# `../data/derived/popPars_100.rds`

# load population params
popPars <- readRDS("../data/derived/popPars_100.rds")

# add IIV on ka and CL/VmaxH
set.seed(192898)
iVmaxH <- rlnorm(100, meanlog=log(40), sdlog=0.2)
ika <- rlnorm(100, meanlog=log(0.849), sdlog=0.2)

# add to popPars
popPars2 <- lapply(1:length(popPars), function(i){
  pars <- c(popPars[[i]],list(VmaxH = iVmaxH[i]))
  pars <- c(pars, list(ka = ika[i]))
  return(pars)
})

# simulate
modA2 <- mread("../model/voriPBPK2")
Kpmu <- 0.56
modA2 <- param(modA2, Kpmu=Kpmu)

#set simulation conditions
bw   <- 73
amt  <- 4*bw
rate <- 4*bw
cmt  <- "VEN"
ii   <- 12
addl <- 13
ss   <- 1
delta <- 0.1
end <- 12

# prepare simulation dataset
idata <- popPars2 %>% bind_rows()

e <- ev(amt=amt, cmt=cmt, ii=ii, addl=addl, rate=rate, ss=ss)

data <- e %>% 
  as_tibble %>% 
  bind_cols(idata) %>%
  mutate(amt = 4*BW, rate = 4*BW) %>%
  select(ID, everything())

#run simulation
system.time(sims <- modA2 %>% 
              mrgsim_d(data=data, delta = delta, end = end, obsonly=T, outvars = c("CP")) %>% 
              filter(row_number() != 1))

# get summary stats for population prediction
hi95 <- function(x) quantile(x, probs = c(0.95))
lo05 <- function(x) quantile(x, probs = c(0.05))

sims2 <- sims %>%
  group_by(time) %>%
  mutate(loCP = lo05(CP),
         medCP = median(CP),
         hiCP = hi95(CP)) %>%
  ungroup() %>%
  filter(ID == first(ID))

#plot population predictions 
gp <- ggplot(data = sims2, aes(x=time)) + 
  geom_line(aes(y=medCP), col="black") +
  geom_ribbon(aes(ymin=loCP, ymax=hiCP), alpha = 0.5) +
  scale_y_continuous(trans = "log10") +
  labs(title="Population simulation", x="time (h)", y="Plasma concentration (mg/L)") +
  theme_bw()
gp

################################################################################
################################################################################

################################################################################
################################## Chunk 16 ####################################
################################################################################

# Use `PKNCA` package to calculate summary statistics for the population

library(PKNCA)

# create dose object
dose_obj <-
  PKNCAdose(
    data,
    amt~time|ID
  )

# create conc object
conc_obj <-
  PKNCAconc(
    sims,
    CP~time|ID
  )

# combine
data_obj <- PKNCAdata(conc_obj, dose_obj)

## Calculate the NCA parameters
results_obj <- pk.nca(data_obj)

## Summarize the results
summary(results_obj)


################################################################################
################################################################################

################################################################################
################################## Chunk 17 ####################################
################################################################################

# Compare doses 3 and 4 mg/kg for adults. Which dose keeps voriconazole concentration 
# at steady state above MIC (1 mg/L) for >= 90% of subjects?

library(cowplot)

data_3mg <- e %>% 
  as_tibble %>% 
  bind_cols(idata) %>%
  mutate(amt = 3*BW, rate = 3*BW) %>%
  select(ID, everything())

#run simulation
sims_3mg <- modA2 %>% 
  mrgsim_d(data=data_3mg, delta = delta, end = end, obsonly=T, outvars = c("CP")) %>% 
  filter(row_number() != 1)

sims2_3mg <- sims_3mg %>%
  group_by(time) %>%
  mutate(loCP = lo05(CP),
         medCP = median(CP),
         hiCP = hi95(CP)) %>%
  ungroup() %>%
  filter(ID == first(ID))

#plot population predictions 
gp_3mg <- ggplot(data = sims2_3mg, aes(x=time)) + 
  geom_line(aes(y=medCP), col="black") +
  geom_ribbon(aes(ymin=loCP, ymax=hiCP), alpha = 0.5) +
  scale_y_continuous(trans = "log10", limits=c(0.3,7)) +
  labs(title="Population simulation (3 mg/kg)", x="time (h)", y="Plasma concentration (mg/L)") +
  theme_bw() + 
  geom_hline(aes(yintercept = 1), lty=2) 

gp_4mg <- gp + 
  labs(title="Population simulation (4 mg/kg)", x="time (h)", y="Plasma concentration (mg/L)") +
  geom_hline(aes(yintercept = 1), lty=2) +
  scale_y_continuous(trans = "log10", limits=c(0.3,7))

gp_34 <- plot_grid(gp_3mg, gp_4mg, ncol=2) 
gp_34

# calculate stats
stats_4mg <- sims %>%
  filter(time != 0) %>%
  group_by(ID) %>%
  mutate(Cmin = min(CP)) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(FLG = ifelse(Cmin < 1, 1, 0))
percLow_4mg <- (nrow(stats_4mg[stats_4mg$FLG==1,]) / nrow(stats_4mg)) * 100

stats_3mg <- sims_3mg %>%
  filter(time != 0) %>%
  group_by(ID) %>%
  mutate(Cmin = min(CP)) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(FLG = ifelse(Cmin < 1, 1, 0))
percLow_3mg <- (nrow(stats_3mg[stats_3mg$FLG==1,]) / nrow(stats_3mg)) * 100 

stats_df <- tibble(Dose = c("3 mg/kg","4 mg/kg"),
                   `Percent of subjects with SS Cmin < MIC` = c(percLow_3mg, percLow_4mg))

stats_df %>%
  knitr::kable() %>%
  kable_styling()

################################################################################
################################################################################

################################################################################
################################## Chunk 18 ####################################
################################################################################

## parallelization ##

# Use the package `mrgsim.apply` (vignette: https://github.com/kylebaron/mrgsim.parallel) 
# to run the same population simulation over parallel cores. How long does it take 
# compared to the single core simulation?

library(parallel)
library(future)
library(mrgsim.parallel)

# set parallelization options
options(future.fork.enable=TRUE, mc.cores = 6L)
plan(multiprocess, workers = 6L)

#run simulation
system.time(sims <- modA2 %>% 
              future_mrgsim_d(data=data, nchunk = 6L, delta = delta, end = end, obsonly=T, outvars = c("CP")) %>% 
              filter(row_number() != 1))

# get summary stats for population
hi95 <- function(x) quantile(x, probs = c(0.95))
lo05 <- function(x) quantile(x, probs = c(0.05))

sims2 <- sims %>%
  group_by(time) %>%
  mutate(loCP = lo05(CP),
         medCP = median(CP),
         hiCP = hi95(CP)) %>%
  ungroup() %>%
  filter(ID == first(ID))

#plot population predictions 
gp <- ggplot(data = sims2, aes(x=time)) + 
  geom_line(aes(y=medCP), col="black") +
  geom_ribbon(aes(ymin=loCP, ymax=hiCP), alpha = 0.5) +
  scale_y_continuous(trans = "log10") +
  labs(title="Population simulation - parallel", x="time (h)", y="Plasma concentration (mg/L)") +
  theme_bw()
gp

################################################################################
################################################################################

#------------------------------------------------------------------------------#

### Shiny interactive document ###

################################################################################
################################## Chunk 19 ####################################
################################################################################

# Explore the simple shiny app saved as `app.R`

library(shiny)
runApp("app.R")

################################################################################
################################################################################


