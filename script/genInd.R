#This function generates the desired individual parameters; using linear interpolaton (PK-Sim)
#source is mainly the PK-Sim Github page https://github.com/Open-Systems-Pharmacology/OSPSuite.Documentation/wiki/Create-Individual-Algorithm
library(tidyverse)
library(Runuran)
library(nloptr)

##load databases
nhanesData <- read.csv("../data/source/dataset2.csv")  #load nhanes anthropometric data
icrpData <- read.csv("../data/source/icrpParams2.csv")  #load physiology data from ICRP publication 89, no separation between gut wall and lumen
SF <- read.csv("../data/source/allometricSF2.csv")  #allometric scaling factors; source: https://github.com/Open-Systems-Pharmacology/OSPSuite.Documentation/wiki/Create-Individual-Algorithm
normSD <- read.csv("../data/source/normalOrgansData2.csv")  #normally distributed organs' sd
lnormSD <- read.csv("../data/source/lnormalOrgansData.csv")  #log normally distributed organs' geometric sd
flow <- read.csv("../data/source/organBloodFlow2.csv")  #organ blood flows
BC <- read.csv("../data/source/organBloodCont2.csv")  #organ blood content

genInd <- function(age, is.male, bw_targ, ht_targ, optimize=TRUE){
  #nhanesData is nhanes anthropometric dataset; icrpData are physiological parameters from ICRP; SF is allomeric scaling factors
  #BC is organ relative blood content
  
  bmi_targ <- bw_targ/ht_targ^2
  
  ##filter nhanes database for desired age and gender
  if(is.male){
    sex <- 1
    dat <- nhanesData %>% dplyr::filter(AGE_YR==age, SEX==sex)  #filter male nhanes data
    df <- icrpData %>% dplyr::select(grep("_m", names(icrpData)))  #filter male physiological data
    names(df) <- gsub("_m", "", names(df))  #remove "_m" from df names
    normSD[,"cv_f"] <- NULL  #remove female CV column
    names(normSD) <- gsub("_m", "", names(normSD))  #remove "_m" from df names
    flow[,"flowPerc_f"] <- NULL  #remove female flows column
    names(flow) <- gsub("_m", "", names(flow))  #remove "_m" from df names
    BC[,"bloodPerc_f"] <- NULL  #remove female blood content column
    names(BC) <- gsub("_m", "", names(BC))  #remove "_m" from df names
  }else{
    sex <- 2
    dat <- nhanesData %>% dplyr::filter(AGE_YR==age, SEX==sex)  #filter female nhanes data
    df <- icrpData %>% dplyr::select(grep("_f", names(icrpData)))  #filter female physiological data
    names(df) <- gsub("_f", "", names(df))  #remove "_f" from df names
    normSD[,"cv_m"] <- NULL  #remove female CV column
    names(normSD) <- gsub("_f", "", names(normSD))  #remove "_m" from df names
    flow[,"flowPerc_m"] <- NULL  #remove male flows column
    names(flow) <- gsub("_f", "", names(flow))  #remove "_m" from df names
    BC[,"bloodPerc_m"] <- NULL  #remove male blood content column
    names(BC) <- gsub("_f", "", names(BC))  #remove "_f" from df names
  }
  
  #get ranges for body weight and height correlated with age and sex
  rangeBW <- quantile(dat$BW, c(0.025, 0.975))  #get the 95% range of body weights
  rangeHT <- quantile(dat$HT/100, c(0.025, 0.975))  #get the 95% range of heights
  rangeBMI <- quantile(dat$BMI, c(0.025, 0.975))  #get the 95% range of body mass index
  
  #throw error when target measurements are out of range
  if (bw_targ < rangeBW[1] || bw_targ > rangeBW[2] || ht_targ < rangeHT[1] || ht_targ > rangeHT[2] || bmi_targ < rangeBMI[1] || bmi_targ > rangeBMI[2]){
    stop("Input values are out of range")
  }
  
  ##get mean bw, ht and bmi
  bw_mean <- exp(mean(log(dat$BW)))  #geomteric mean for lognormally distributed weights
  ht_mean <- mean(dat$HT)/100  #arithmetic mean for normally distributed heights in m
  bmi_mean <- exp(mean(log(dat$BMI)))  #geomteric mean for lognormally distributed BMI; Jim and I think this is better   
  
  ##get mean physiological parameters scaled by linear interpolation with height
  df2 <- df %>% dplyr::select(-c(bw, ht, bsa))  #get rid of anthropometric measurements
  ht_ref <- df$ht/100  #get reference body heights from ICRP in m
  linIntFns <- apply(df2, 2, FUN=function(x) approxfun(ht_ref, x, rule=2))  #linear interpolation functions with height; if outside range use closest value
  physPars <- lapply(linIntFns, FUN=function(x) x(ht_mean))  #get list of mean physiological parameters
  
  ##########piling mean physilogical values and standard deviations############
  df_temp <- data.frame(organ=names(physPars), means=as.numeric(physPars))  #get scaling for blood and cardiac output
  blVol <- df_temp$means[df_temp$organ == "bl"]  #extract blood volume
  df_temp <- bind_rows(df_temp, data.frame(organ=c("ve","ar"), means=c(0, 0)))  #temporarily use whole blood volume for arterial and venous blood vols
  df_temp <- df_temp %>% dplyr::filter(organ != "bl")  #remove blood volume from df
  l_temp <- split(df_temp, df_temp$organ %in% c("co"))  #separate vols from co
  df_vols <- l_temp[["FALSE"]]  
  df_co <- l_temp[["TRUE"]]  
  
  #################################getting organ volumes##########################################
  #add blood content to vols
  df_vols <- left_join(df_vols, BC, by="organ") 
  df_vols <- df_vols %>% dplyr::mutate(means=means + (bloodPerc*blVol/100))
  df_temp <- dplyr::bind_rows(df_vols %>% dplyr::select(organ, means), df_co)  #join with co again
  
  normOrgan <- normSD$organ
  lnormOrgan <- lnormSD$organ
  df_norm <- df_temp %>% dplyr::filter(organ %in% normOrgan)  #prepare a df for normally distributed organs
  df_norm <- left_join(df_norm, normSD, by="organ")  #integrate with stds
  df_norm <- df_norm %>% dplyr::mutate(std = cv*means/100) %>% dplyr::select(-cv)  #get standard deviations from CVs then drop it
  df_lnorm <- df_temp %>% dplyr::filter(organ %in% lnormOrgan)  #prepare a df for lognormally distributed organs
  df_lnorm <- left_join(df_lnorm, lnormSD, by="organ")  #integrate with stds
  names(df_lnorm)[names(df_lnorm) == "geomSD"] <- "std"  #change geomSD to std
  df3 <- bind_rows(df_norm, df_lnorm)  #get combined dataframe for normally and lognormally distributed organs
  
  ##scaling of physiological parameters
  ht_rel <- ht_targ/ht_mean  #get relative height ratio
  SF2 <- bind_rows(SF, data.frame(organ="co", factor=0.75))  #add allometric scaling factor for cardiac output
  df3 <- left_join(df3, SF2, by=("organ"))  #integrate allometric scaling factors
  df3 <- df3 %>% dplyr::mutate(means_scaled = ifelse(organ %in% normOrgan, means*ht_rel^factor, means + log(ht_rel^factor)))  #get scaled mean values
  df3 <- df3 %>% dplyr::mutate(std_scaled = ifelse(organ %in% normOrgan, std*ht_rel^factor, std))  #get scaled standard deviations
  df4 <- df3 %>% dplyr::filter(organ != "co")  #remove cardiac output
  
  ######################### optimization of organ volumes #################################
  #set initial values for volumes
  l_temp <- split(df4, df4$organ %in% c("ad"))  #separate vols from ad
  df_opt <- l_temp[["FALSE"]]  
  df_ad <- l_temp[["TRUE"]]
  
  #df_opt <- df4 %>% dplyr::filter(organ != "ad") #remove adipose as we will not optimize for it
  pert <- runif(1)  #perturbation to parameters
  df_opt <- df_opt %>% dplyr::mutate(pertMeans = ifelse(organ %in% normOrgan,
                                                 means_scaled + std_scaled*pert,
                                                 means_scaled*exp(log(std_scaled)*pert))) %>%
                       dplyr::mutate(pertMeans = ifelse(organ == "sk", pertMeans*sqrt((bw_targ/(bw_mean*ht_rel^2))), pertMeans))#temporary fix
  #exp(log(means_scaled^2/(sqrt(std_scaled^2 + means_scaled^2))) + sqrt(log(std_scaled^2/(means_scaled^2) + 1))*pert))) #get perturbed volumes ##formulae from https://en.wikipedia.org/wiki/Log-normal_distribution
  ad <- bw_targ - sum(df_opt$pertMeans) #compute adipose volume by subtracting current weight from target bw_targ
  p_ad <- plnorm(ad, meanlog=log(df_ad$means_scaled), sdlog=log(df_ad$std_scaled))  #compute prob for adipose vol
  
  optVols <- df_opt$pertMeans  #initial parameter values
  
  #set the optimization function to return probability
  optimVols <- function(optVols){
    df_opt$pertMeans <- optVols
    ad <- bw_targ - sum(df_opt$pertMeans)
    p_ad <- dlnorm(ad, meanlog=log(df_ad$means_scaled), sdlog=log(df_ad$std_scaled))
    skFactor <- sqrt(bw_targ/(bw_mean*ht_rel^2)) 
    df_opt <- df_opt %>% dplyr::mutate(means_scaled = ifelse(organ %in% c("sk"), 
                                                   means_scaled + skFactor, 
                                                   means_scaled)) #get a new distribution for skin
    df_opt <- df_opt %>% dplyr::mutate(probs = ifelse(organ %in% normOrgan, dnorm(optVols, mean=means_scaled, sd=std_scaled),
                                               dlnorm(optVols, meanlog=log(means_scaled), sdlog=log(std_scaled)))) #needs to change
    p <- -sum(c(log(df_opt$probs), log(p_ad)))  #optimize for the negative log probability because we want to maximize not minimize
    #p <- -prod(c(df_opt$probs, p_ad))
    return(p)
  }
  
  #optMod <- optim(optVols, optimVols, method="Nelder-Mead")  #optimize probabilities; optim faster than nloptr neldermead
  #optMod <- neldermead(optVols, optimVols)  #Nelder-Mead simplex; nloptr package
  if(is.optimized){
    optMod <- nloptr::newuoa(optVols, optimVols)
    df_opt <- df_opt %>% dplyr::mutate(optimized=optMod$par)
  }else{
    df_opt <- df_opt %>% dplyr::mutate(optimized=means_scaled)
  }
  
  #ad <- bw_targ - sum(df_opt$optimized)  #still need to be corrected 

  l_ov <- as.list(df_opt$optimized)  #get the final scaled organ volumes in a list
  names(l_ov) <- paste("V", df_opt$organ, sep="")  #name the list
  l_ov <- c(l_ov, Vad=ad)
  
  #################################getting organ blood flows#####################################
  co <- df3$means_scaled[df3$organ == "co"]
  flow2 <- flow %>% dplyr::mutate(bfs = (flowPerc*co/100)*60) %>% dplyr::filter(!organ %in% c("ve","ar","lu","li")) #get flow rates in L/h and remove lungs, liver, arterial and venous blood flows as they will be calculated within the model
  a <- ifelse(runif(1) > 0.5, 0, -Inf)
  b <- ifelse(a == 0, Inf, 0)
  flow2 <- flow2 %>% dplyr::mutate(bfs = bfs + urnorm(length(bfs), lb=a, ub=b, mean=0, sd=0.05*abs(bfs)))
  l_bf <- as.list(flow2$bfs)  #list of blood flows
  names(l_bf) <- paste("Q", flow2$organ, sep="")  #name the list
  
  ##getting final parameter list
  pars <- c(l_ov, l_bf, BW=bw_targ, HT=ht_targ, BMI=bmi_targ, SEX=sex)
}


