## This R script compiles and runs simulations for voriPBPK model; tasks here follow `hands-on_voriPBPK`

## load libraries
rm(list=ls())
library(dplyr)
library(ggplot2)
library(mrgsolve)


## adjust general theme for plotting
th <- theme(plot.title=element_text(size=30),
            axis.title=element_text(size=25),
            axis.text=element_text(size=20),
            legend.text=element_text(size=20))

#############################################################################################
##############################  Chunk 1: Compile model  #####################################
#############################################################################################

modA <- mread("voriPBPK", "models")

#############################################################################################
#############################################################################################


#############################################################################################
#########################  Chunk 2: Run a simple simulation  ################################
#############################################################################################

## Simulate a 100 mg dose given as an IV bolus dose
modA %>%
  ev(amt=100, cmt="VEN") %>%
  mrgsim() %>%
  plot()

#############################################################################################
#############################################################################################


#############################################################################################
#################################  Chunk 3: Task 2  #########################################
#############################################################################################
## Simulate the following scenarios:
## A) 4 mg/kg voriconazole IV infusion dosing in an adult male (73 kg) infused over an hour 
## twice a day for 7 days

### Code here:


## B) 200 mg voriconazole given orally to adult male twice a day for 7 days

### Code here:


#############################################################################################
#############################################################################################


#############################################################################################
#################################  Chunk 4: Task 3  #########################################
#############################################################################################
## Compare the steady state plasma drug concentration-time profiles from previous simulations 
## to the observed data

## A) IV infusion (observed data in `inst/data/Adult_IV.csv`)
obs <- read.csv("inst/data/Adult_IV.csv")  #load observed data

### Run this ###
wt <- 73  #adult body weight
dose <- 4*wt  
rate <- 4*wt
cmt <- "VEN"  #intravenous infusion

# simulate
sim <- as.data.frame(modA %>% 
                       ev(amt=dose, cmt=cmt, ii=12, addl=13, rate=rate, ss=1) %>% 
                       mrgsim(delta = 0.1, end = 12)) %>% 
  dplyr::filter(row_number() != 1)  

# plot
gp <- ggplot() + 
  geom_point(data = obs, aes(x=time, y=obs, col="observed"), size=2.5) + 
  geom_errorbar(data = obs, aes(x = time, y = obs, ymin=obs-sd, ymax=obs+sd), width=0) +
  geom_line(data = sim, aes(x=time, y=Cvenous, col="sim"), lwd=1) + 
  scale_colour_manual(name='', 
                      values=c('sim'='black', 'observed'='black'), 
                      breaks=c("observed","sim"),
                      labels=c("observed","predicted")) +
  guides(colour = guide_legend(override.aes = list(linetype=c(0,1), shape=c(16, NA)))) +
  labs(x="time (h)", y="Plasma concentration (mg/L)") +
  th
gp


## B) Oral dosing (observed data in `inst/data/Adult_PO.csv`)
obs <- read.csv("inst/data/Adult_PO.csv")  #load observed data

wt <- 73
dose <- 200  
cmt <- "GUTLUMEN"  #intravenous infusion

# simulate
sim <- as.data.frame(modA %>%
                       ev(amt=dose, cmt=cmt, ii=12, addl=13, ss=1) %>% 
                       mrgsim(delta = 0.1, end = 12)) %>% 
  dplyr::filter(row_number() != 1)  

# plot
gp <- ggplot() + 
  geom_point(data = obs, aes(x=time, y=obs, col="observed"), size=2.5) + 
  geom_errorbar(data = obs, aes(x = time, y = obs, ymin=obs-sd, ymax=obs+sd), width=0) +
  geom_line(data = sim, aes(x=time, y=Cvenous, col="sim"), lwd=1) + 
  scale_colour_manual(name='', 
                      values=c('sim'='black', 'observed'='black'), 
                      breaks=c("observed","sim"),
                      labels=c("observed","predicted")) +
  guides(colour = guide_legend(override.aes = list(linetype=c(0,1), shape=c(16, NA)))) +
  labs(x="time (h)", y="Plasma concentration (mg/L)") +
  th
gp

#############################################################################################
#############################################################################################

#############################################################################################
#################################  Chunk 5: Task 5  #########################################
#############################################################################################
## Compile the `voriPBPK_absorption` model and re-run the previous oral simulation 
modA_abs <- mread("voriPBPK_absorption", "models")  #compile model

obs <- read.csv("inst/data/Adult_PO.csv")  #load observed data

wt <- 73
dose <- 200  
cmt <- "GUTLUMEN"  #intravenous infusion

# simulate
sim <- as.data.frame(modA_abs %>%
                       ev(amt=dose, cmt=cmt, ii=12, addl=13, ss=1) %>% 
                       mrgsim(delta = 0.1, end = 12)) %>% 
  dplyr::filter(row_number() != 1)  

# plot
gp <- ggplot() + 
  geom_point(data = obs, aes(x=time, y=obs, col="observed"), size=2.5) + 
  geom_errorbar(data = obs, aes(x = time, y = obs, ymin=obs-sd, ymax=obs+sd), width=0) +
  geom_line(data = sim, aes(x=time, y=Cvenous, col="sim"), lwd=1) + 
  scale_colour_manual(name='', 
                      values=c('sim'='black', 'observed'='black'), 
                      breaks=c("observed","sim"),
                      labels=c("observed","predicted")) +
  guides(colour = guide_legend(override.aes = list(linetype=c(0,1), shape=c(16, NA)))) +
  labs(x="time (h)", y="Plasma concentration (mg/L)") +
  th
gp

#############################################################################################
#############################################################################################

#############################################################################################
#################################  Chunk 6: Task 6  #########################################
#############################################################################################
## Run simple graphical sensitivity analyses for absorption parameters: permeability, 
## intestinal transit time and solubility

##' Define an intervention
e <- ev(amt = 200, cmt = "GUTLUMEN", ii = 12, addl = 13, ss = 1)

##' Sensitivity analysis on drug solubility (S_lumen)
idata <- data_frame(S_lumen = c(0.39*500, 0.39*1000, 0.39*2000)) %>% 
  mutate(ID = seq(n()))

modA_abs %>% 
  carry_out(S_lumen) %>% 
  mrgsim_ei(e,idata, delta = 0.1, end = 12, recsort=3, obsonly=TRUE) %>%
  plot(Cvenous ~ time, 
       group = S_lumen, 
       auto.key = TRUE,
       main=list(label=expression(S["int"]~(mg/mL)), cex=1), 
       xlab=list(label="time (h)", cex=0.75),
       ylab=list(label="Plasma concentration (mg/L)", cex=0.75), 
       scale=list(cex=0.5), 
       lty=1:3, 
       col="black",
       #auto.key=list(cex=2, lines=T, points=F),
       key=list(space="top",
                lines=list(col=c("black","black","black"), lty=1:3),
                text=list(c("0.195","0.39","0.78"), cex=0.5)
       ))


##' Sensitivity analysis on permeability (fperm)
idata <- data_frame(fperm=c(0.5, 1, 2)) %>% 
  mutate(ID = seq(n()))

modA_abs %>% 
  carry_out(fperm) %>% 
  mrgsim_ei(e,idata, delta = 0.1, recsort=3, obsonly=TRUE, end = 12) %>%
  plot(Cvenous ~ time, 
       group = fperm, 
       auto.key = TRUE,
       main=list(label=expression(paste("P"[eff], " x10",phantom()^{-4}," (cm/s)")), cex=1),
       xlab=list(label="time (h)", cex=0.75),
       ylab=list(label="Plasma concentration (mg/L)", cex=0.75), 
       scale=list(cex=0.5), 
       lty=1:3, 
       col="black",
       #auto.key=list(cex=2, lines=T, points=F),
       key=list(space="top",
                lines=list(col=c("black","black","black"), lty=1:3),
                text=list(c("0.073","0.145","0.29"), cex=0.5)
       ))


##' Sensitivity analysis on intestinal transit time (ITT)
idata <- data_frame(ITT=c(3.32/2, 3.32, 3.32*2)) %>% 
  mutate(ID = seq(n()))

modA_abs %>% 
  carry_out(ITT) %>% 
  mrgsim_ei(e,idata, delta = 0.1, recsort=3, obsonly=TRUE, end = 12) %>%
  plot(Cvenous ~ time, 
       group = ITT, 
       auto.key = TRUE,
       main=list(label=expression(ITT~(h)), cex=1),
       xlab=list(label="time (h)", cex=0.75),
       ylab=list(label="Plasma concentration (mg/L)", cex=0.75),
       scale=list(cex=0.5),
       lty=1:3,
       col="black",
       #auto.key=list(cex=2, lines=T, points=F),
       key=list(space="top",
                lines=list(col=c("black","black","black"), lty=1:3),
                text=list(c("1.66","3.32","6.64"), cex=0.5)
       ))


##' Sensitivity analysis on intestinal clearance (MPPGI)
idata <- data_frame(MPPGI=c(1.212/2, 1.212, 1.212*2)) %>% 
  mutate(ID = seq(n()))

modA_abs %>% 
  carry_out(MPPGI) %>% 
  mrgsim_ei(e,idata, delta = 0.1, recsort=3, obsonly=TRUE, end = 12) %>%
  plot(Cvenous ~ time, 
       group = MPPGI, 
       auto.key = TRUE,
       main=list(label=expression(Cl["Gu"]~(mL/min/kg)), cex=1),
       xlab=list(label="time (h)", cex=0.75),
       ylab=list(label="Plasma concentration (mg/L)", cex=0.75),
       scale=list(cex=0.5),
       lty=1:3,
       col="black",
       #auto.key=list(cex=2, lines=T, points=F),
       key=list(space="top",
                lines=list(col=c("black","black","black"), lty=1:3),
                text=list(c("0.035","0.07","0.14"), cex=0.5)
       ))

#############################################################################################
#############################################################################################

#############################################################################################
#################################  Chunk 7: Task 7  #########################################
#############################################################################################
## Finetune the influential absorption parameters and re-run oral simulation
## Update model
modA_abs <- param(modA_abs, fperm=0.47, MPPGI=30.3/25)

## load observed data
obs <- read.csv("inst/data/Adult_PO.csv")  #load observed data

wt <- 73
dose <- 200  
cmt <- "GUTLUMEN"  #intravenous infusion

# simulate
sim <- as.data.frame(modA_abs %>%
                       ev(amt=dose, cmt=cmt, ii=12, addl=13, ss=1) %>% 
                       mrgsim(delta = 0.1, end = 12)) %>% 
  dplyr::filter(row_number() != 1)  

# plot
gp <- ggplot() + 
  geom_point(data = obs, aes(x=time, y=obs, col="observed"), size=2.5) + 
  geom_errorbar(data = obs, aes(x = time, y = obs, ymin=obs-sd, ymax=obs+sd), width=0) +
  geom_line(data = sim, aes(x=time, y=Cvenous, col="sim"), lwd=1) + 
  scale_colour_manual(name='', 
                      values=c('sim'='black', 'observed'='black'), 
                      breaks=c("observed","sim"),
                      labels=c("observed","predicted")) +
  guides(colour = guide_legend(override.aes = list(linetype=c(0,1), shape=c(16, NA)))) +
  labs(x="time (h)", y="Plasma concentration (mg/L)") +
  th
gp

#############################################################################################
#############################################################################################

#############################################################################################
#################################  Chunk 8: Task 8  #########################################
#############################################################################################
## Generate pediatric model with the integrated absorption model

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
                Qgu = 0.16*3.4*60, 
                Qhe = 0.04*3.4*60,
                Qki = 0.19*3.4*60,
                Qmu = 0.17*3.4*60,
                Qsp = 0.03*3.4*60,
                Qha = 0.065*3.4*60, 
                Qlu = 3.4*60,
                MPPGL = 26,
                VmaxH = 120.5,
                KmH = 11,
                MPPGI = 26/25,
                VmaxG = 120.5,
                KmG = 11,
                L = 170)

### Code here:


#############################################################################################
#############################################################################################

#############################################################################################
#################################  Chunk 9: Task 9  #########################################
#############################################################################################
## Simulate the following scenarios:
## A) 4 mg/kg voriconazole IV infusion dosing in a male child subject infused 
## with a rate of 3 mg/kg/h twice a day for 7 days

### Code here:


## B) 4 mg/kg given orally twice a day to male child for 7 days

### Code here:


#############################################################################################
#############################################################################################

#############################################################################################
#################################  Chunk 10: Task 10  #######################################
#############################################################################################
## Compare the steady state plasma drug concentration-time profiles from previous simulations 
## to the observed data

## A) IV infusion (observed data in `inst/data/Pediatric_IV.csv`)
obs <- read.csv("inst/data/Pediatric_IV.csv")  #load observed data

### Run this:
wt <- 19  #adult body weight
dose <- 4*wt  
rate <- 3*wt
cmt <- "VEN"  #intravenous infusion

# simulate
sim <- as.data.frame(modP_abs %>% 
                       ev(amt=dose, cmt=cmt, ii=12, addl=13, rate=rate, ss=1) %>% 
                       mrgsim(delta = 0.1, end = 12)) %>% 
  dplyr::filter(row_number() != 1)  

# plot
gp <- ggplot() + 
  geom_point(data = obs, aes(x=time, y=obs, col="observed"), size=2.5) + 
  geom_errorbar(data = obs, aes(x = time, y = obs, ymin=obs-sd, ymax=obs+sd), width=0) +
  geom_line(data = sim, aes(x=time, y=Cvenous, col="sim"), lwd=1) + 
  scale_colour_manual(name='', 
                      values=c('sim'='black', 'observed'='black'), 
                      breaks=c("observed","sim"),
                      labels=c("observed","predicted")) +
  guides(colour = guide_legend(override.aes = list(linetype=c(0,1), shape=c(16, NA)))) +
  labs(x="time (h)", y="Plasma concentration (mg/L)") +
  th
gp


## B) Oral dosing (observed data in `inst/data/Pediatric_PO.csv`)
obs <- read.csv("inst/data/Pediatric_PO.csv")  #load observed data

wt <- 19
dose <- 4*wt  
cmt <- "GUTLUMEN"  #intravenous infusion

# simulate
sim <- as.data.frame(modP_abs %>%
                       ev(amt=dose, cmt=cmt, ii=12, addl=7, ss=1) %>% 
                       mrgsim(delta = 0.1, end = 12)) %>% 
  dplyr::filter(row_number() != 1)  

# plot
gp <- ggplot() + 
  geom_point(data = obs, aes(x=time, y=obs, col="observed"), size=2.5) + 
  geom_errorbar(data = obs, aes(x = time, y = obs, ymin=obs-sd, ymax=obs+sd), width=0) +
  geom_line(data = sim, aes(x=time, y=Cvenous, col="sim"), lwd=1) + 
  scale_colour_manual(name='', 
                      values=c('sim'='black', 'observed'='black'), 
                      breaks=c("observed","sim"),
                      labels=c("observed","predicted")) +
  guides(colour = guide_legend(override.aes = list(linetype=c(0,1), shape=c(16, NA)))) +
  labs(x="time (h)", y="Plasma concentration (mg/L)") +
  th
gp

#############################################################################################
#############################################################################################

#############################################################################################
#################################  Chunk 11: Task 11  #######################################
#############################################################################################
## Run shiny app and play with it

#############################################################################################
#############################################################################################