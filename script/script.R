# load libraries
library(tidyverse)
library(mrgsolve)
library(FME)
library(sensitivity)
library(nloptr)

# Chunks 1-2 are relevant to simplePBPK
# Chunks 3-10 are relevant to voriPBPK model
# Chunks 11-13 are relevant to sensitivity analysis


#------------------------------------------------------------------------------#

### Simple PBPK ###

################################################################################
################################## Chunk 1  ####################################
################################################################################

# Compile simplePBPK model
mod <- mread("simplePBPK", "models")

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

# Run graphical sensitivity analysis for the muscle:plasma (`Kpmu`) and 
# lung:plasma (`Kplu`) partition coefficients using adult IV data.



################################################################################
################################################################################



