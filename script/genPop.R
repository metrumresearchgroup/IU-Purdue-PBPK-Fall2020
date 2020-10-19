#This function generates the desired individual parameters. It takes in the model (mod), number of subjects (n), ranges for age, height 
#and weight and and the percentage of females in the population 
.libPaths("lib")
library(tidyverse)
source("genInd.R")

genPop <- function(nSubj, minAge, maxAge, femPerc, minBW, maxBW, minHT, maxHT){
  nFemale <- (femPerc/100)*nSubj
  nMale <- nSubj-nFemale
  pars_m <- rep(list(), nMale)
  pars_f <- rep(list(), nFemale)
  
  ind <- 1
  while(ind <= nMale){
    age <- round(runif(1, minAge, maxAge))
    is.male <- TRUE
    dat <- nhanesData %>% dplyr::filter(AGE_YR == age & SEX == 1)
    bw_targ <- rlnorm(1, mean=mean(log(dat$BW)), sd=sd(log(dat$BW)))
    ht_targ <- (rnorm(1, mean=mean(dat$HT), sd=sd(dat$HT)))/100
    bmi_targ <- bmi_targ <- bw_targ/ht_targ^2

    t <- try(pars_m[[ind]] <- genInd(age=age, is.male=is.male, bw_targ=bw_targ, ht_targ=ht_targ))  #get the individual parameters
    if("try-error" %in% class(t)){
      ind <- ind
    }else{
      ind <- ind + 1
    }
  }
  
  ind <- 1
  while(ind <= nFemale){
    age <- round(runif(1, minAge, maxAge))
    is.male <- FALSE
    dat <- nhanesData %>% dplyr::filter(AGE_YR == age & SEX == 2)
    bw_targ <- rlnorm(1, mean=mean(log(dat$BW)), sd=sd(log(dat$BW)))
    ht_targ <- (rnorm(1, mean=mean(dat$HT), sd=sd(dat$HT)))/100
    bmi_targ <- bmi_targ <- bw_targ/ht_targ^2
    
    t <- try(pars_f[[ind]] <- genInd(age=age, is.male=is.male, bw_targ=bw_targ, ht_targ=ht_targ))  #get the individual parameters and get error if input out of range and try again
    if("try-error" %in% class(t)){
      ind <- ind
    }else{
      ind <- ind + 1
    }
  }
  
  pars <- c(pars_m, pars_f)
  return(pars)
}

