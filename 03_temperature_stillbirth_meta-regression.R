##############################################################################
#
# Contact: Jeroen de Bont (jeroen.de.bont@ki.se)
#
##############################################################################
#
# Script:  Meta-regression between temperature (lag-06) and stillbirth
#
##############################################################################

#rm(list=ls())
library(readxl)
library(foreign)
library(data.table)
library(lubridate)
library(gdata)
library(rSPARCS)
library(dlnm)
library(dplyr)
library(survival)
library(splines)
library(mixmeta)

#########################################################################################################
##### Prepare the dataset
#########################################################################################################

# setnames
#setwd("/Users/jeroen.de.bont/Desktop/3_paper_alert/")
setwd("/Volumes/projects$/EnvEpi-Projekt/CHAMNHA/4_paper_alert")

##### Open the final analyses data to select the IDs
alert <- as.data.table(readRDS("03_data/alert_cc.rds"))

# check missing
cc2 <-alert

cc2 <- cc2[order(alert_id, day)]

cc2[is.na(temp_s0)] # no missing
cc2[is.na(intra_still)] # more missing as if you have ante_Still birth you remove them
cc2[is.na(ante_still)] # no missing
cc2[is.na(per_mort)] # no missing

#####################################################

# we are using the 9x9km resolution models which is labelled temp_s0_9. We changed the labelling
cc2[, c("temp_s0","temp_s1","temp_s2","temp_s3","temp_s4","temp_s5","temp_s6") := NULL]

setnames(cc2, c("temp_s0_9","temp_s1_9","temp_s2_9","temp_s3_9","temp_s4_9","temp_s5_9","temp_s6_9"),
         c("temp_s0","temp_s1","temp_s2","temp_s3","temp_s4","temp_s5","temp_s6"))

# keep the same variables as the previous version
cc2 <- as.data.frame(cc2)
cc2$hosp <- as.numeric(cc2$hosp)
cc.still.ben <- subset(cc2, all_still == "Yes" & country == "Benin")
cc.intra.ben <- subset(cc2, intra_still == "Yes" & country == "Benin")
cc.ante.ben <- subset(cc2, ante_still == "Yes" & country == "Benin")
cc.per.ben <- subset(cc2, per_mort == "Yes" & country == "Benin")

cc.still.mal <- subset(cc2, all_still == "Yes" & country == "Malawi")
cc.intra.mal <- subset(cc2, intra_still == "Yes" & country == "Malawi")
cc.ante.mal <- subset(cc2, ante_still == "Yes" & country == "Malawi")
cc.per.mal <- subset(cc2, per_mort == "Yes" & country == "Malawi")

cc.still.tan <- subset(cc2, all_still == "Yes" & country == "Tanzania")
cc.intra.tan <- subset(cc2, intra_still == "Yes" & country == "Tanzania")
cc.ante.tan <- subset(cc2, ante_still == "Yes" & country == "Tanzania")
cc.per.tan <- subset(cc2, per_mort == "Yes" & country == "Tanzania")

cc.still.ug <- subset(cc2, all_still == "Yes" & country == "Uganda")
cc.intra.ug <- subset(cc2, intra_still == "Yes" & country == "Uganda")
cc.ante.ug <- subset(cc2, ante_still == "Yes" & country == "Uganda")
cc.per.ug <- subset(cc2, per_mort == "Yes" & country == "Uganda")


##-------------------------------------------------------------------------------------####
## 1. FUNCTION TO ESTIMATE MINIMUM OF A EXPOSURE-RESPONSE FUNCTION FROM A FITTED MODEL ####
##-------------------------------------------------------------------------------------####

findmin <- function(basis, model=NULL, coef=NULL, vcov=NULL, at=NULL, from=NULL,
                    to=NULL, by=NULL, sim=FALSE, nsim=5000) {
  
  ################################################################################
  # R code from https://github.com/gasparrini/2017_tobias_Epidem_Rcodedata/blob/master/findmin.R
  # 
  #   ARGUMENTS:
  #   - basis: A SPLINE OR OTHER BASIS FOR AN EXPOSURE x CREATED BY DLNM FUNCTION 
  #            CROSSBASIS OR ONEBASIS
  #   - model: THE FITTED MODEL
  #   - coef AND vcov: COEF AND VCOV FOR basis IF model IS NOT PROVIDED
  #
  #   - at: A NUMERIC VECTOR OF x VALUES OVER WHICH THE MINIMUM IS SOUGHT
  #   OR 
  #   - from, to: RANGE OF x VALUES OVER WHICH THE MINIMUM IS SOUGHT.
  #   - by: INCREMENT OF THE SEQUENCES x VALUES OVER WHICH THE MINIMUM IS SOUGHT
  # 
  #   - sim: IF BOOTSTRAP SIMULATION SAMPLES SHOULD BE RETURNED
  #   - nsim: NUMBER OF SIMULATION SAMPLES
  ################################################################################
  
  
  ################################################################################
  # CREATE THE BASIS AND EXTRACT COEF-VCOV
  #
  # CHECK AND DEFINE BASIS  
  if(!any(class(basis)%in%c("crossbasis","onebasis")))
    stop("the first argument must be an object of class 'crossbasis' or 'onebasis'")
  #
  # INFO
  one <- any(class(basis)%in%c("onebasis"))
  attr <- attributes(basis)
  range <- attr(basis,"range")
  if(is.null(by)) by <- 0.1
  lag <- if(one) c(0,0) else cb=attr(basis,"lag")
  if(is.null(model)&&(is.null(coef)||is.null(vcov)))
    stop("At least 'model' or 'coef'-'vcov' must be provided")
  name <- deparse(substitute(basis))
  cond <- if(one) paste(name,"[[:print:]]*b[0-9]{1,2}",sep="") else 
    paste(name,"[[:print:]]*v[0-9]{1,2}\\.l[0-9]{1,2}",sep="")
  #
  # SET COEF, VCOV CLASS AND LINK
  if(!is.null(model)) {
    model.class <- class(model)
    coef <- dlnm:::getcoef(model,model.class)
    ind <- grep(cond,names(coef))
    coef <- coef[ind]
    vcov <- dlnm:::getvcov(model,model.class)[ind,ind,drop=FALSE]
    model.link <- dlnm:::getlink(model,model.class)
  } else model.class <- NA
  #
  # CHECK
  if(length(coef)!=ncol(basis) || length(coef)!=dim(vcov)[1] ||
     any(is.na(coef)) || any(is.na(vcov)))
    stop("model or coef/vcov not consistent with basis")
  #
  # DEFINE at
  at <- dlnm:::mkat(at,from,to,by,range,lag,bylag=1)
  predvar <- if(is.matrix(at)) rownames(at) else at
  predlag <- dlnm:::seqlag(lag,by=1)
  #
  # CREATE THE MATRIX OF TRANSFORMED CENTRED VARIABLES (DEPENDENT ON TYPE)
  type <- if(one) "one" else "cb"
  Xpred <- dlnm:::mkXpred(type,basis,at,predvar,predlag,cen=NULL)
  Xpredall <- 0
  for(i in seq(length(predlag))) {
    ind <- seq(length(predvar))+length(predvar)*(i-1)
    Xpredall <- Xpredall + Xpred[ind,,drop=FALSE]
  }
  #  
  ################################################################################
  # FIND THE MINIMUM
  #
  pred <- drop(Xpredall%*%coef)
  ind <- which.min(pred)
  min <- predvar[ind]
  #
  ################################################################################
  # SIMULATIONS
  #
  if(sim) {
    # SIMULATE COEFFICIENTS
    k <- length(coef)
    eigen <- eigen(vcov)
    X <- matrix(rnorm(length(coef)*nsim),nsim)
    coefsim <- coef + eigen$vectors %*% diag(sqrt(eigen$values),k) %*% t(X)
    # COMPUTE MINIMUM
    minsim <- apply(coefsim,2,function(coefi) {
      pred <- drop(Xpredall%*%coefi)
      ind <- which.min(pred)
      return(predvar[ind])
    })
  }
  #
  ################################################################################
  #
  res <- if(sim) minsim else min
  #
  return(res)
}


##-------------------------------------------------------------------------------------####
## 2. FUNCTION TO ESTIMATE MAXIMUM OF A EXPOSURE-RESPONSE FUNCTION FROM A FITTED MODEL ####
##-------------------------------------------------------------------------------------####

findmax <- function(basis,model=NULL,coef=NULL,vcov=NULL,at=NULL,from=NULL,
                    to=NULL,by=NULL,sim=FALSE,nsim=5000) {
  #
  ################################################################################
  #   Adapted from R code findmin()
  #   ARGUMENTS:
  #   - basis: A SPLINE OR OTHER BASIS FOR AN EXPOSURE x CREATED BY DLNM FUNCTION 
  #            CROSSBASIS OR ONEBASIS
  #   - model: THE FITTED MODEL
  #   - coef AND vcov: COEF AND VCOV FOR basis IF model IS NOT PROVIDED
  #
  #   - at: A NUMERIC VECTOR OF x VALUES OVER WHICH THE MINIMUM IS SOUGHT
  #   OR 
  #   - from, to: RANGE OF x VALUES OVER WHICH THE MINIMUM IS SOUGHT.
  #   - by: INCREMENT OF THE SEQUENCES x VALUES OVER WHICH THE MINIMUM IS SOUGHT
  # 
  #   - sim: IF BOOTSTRAP SIMULATION SAMPLES SHOULD BE RETURNED
  #   - nsim: NUMBER OF SIMULATION SAMPLES
  ################################################################################
  
  
  ################################################################################
  # CREATE THE BASIS AND EXTRACT COEF-VCOV
  #
  # CHECK AND DEFINE BASIS  
  if(!any(class(basis)%in%c("crossbasis","onebasis")))
    stop("the first argument must be an object of class 'crossbasis' or 'onebasis'")
  #
  # INFO
  one <- any(class(basis)%in%c("onebasis"))
  attr <- attributes(basis)
  range <- attr(basis,"range")
  if(is.null(by)) by <- 0.1
  lag <- if(one) c(0,0) else cb=attr(basis,"lag")
  if(is.null(model)&&(is.null(coef)||is.null(vcov)))
    stop("At least 'model' or 'coef'-'vcov' must be provided")
  name <- deparse(substitute(basis))
  cond <- if(one) paste(name,"[[:print:]]*b[0-9]{1,2}",sep="") else 
    paste(name,"[[:print:]]*v[0-9]{1,2}\\.l[0-9]{1,2}",sep="")
  #
  # SET COEF, VCOV CLASS AND LINK
  if(!is.null(model)) {
    model.class <- class(model)
    coef <- dlnm:::getcoef(model,model.class)
    ind <- grep(cond,names(coef))
    coef <- coef[ind]
    vcov <- dlnm:::getvcov(model,model.class)[ind,ind,drop=FALSE]
    model.link <- dlnm:::getlink(model,model.class)
  } else model.class <- NA
  #
  # CHECK
  if(length(coef)!=ncol(basis) || length(coef)!=dim(vcov)[1] ||
     any(is.na(coef)) || any(is.na(vcov)))
    stop("model or coef/vcov not consistent with basis")
  #
  # DEFINE at
  at <- dlnm:::mkat(at,from,to,by,range,lag,bylag=1)
  predvar <- if(is.matrix(at)) rownames(at) else at
  predlag <- dlnm:::seqlag(lag,by=1)
  #
  # CREATE THE MATRIX OF TRANSFORMED CENTRED VARIABLES (DEPENDENT ON TYPE)
  type <- if(one) "one" else "cb"
  Xpred <- dlnm:::mkXpred(type,basis,at,predvar,predlag,cen=NULL)
  Xpredall <- 0
  for(i in seq(length(predlag))) {
    ind <- seq(length(predvar))+length(predvar)*(i-1)
    Xpredall <- Xpredall + Xpred[ind,,drop=FALSE]
  }
  #  
  ################################################################################
  # FIND THE MINIMUM
  #
  pred <- drop(Xpredall%*%coef)
  ind <- which.max(pred)
  min <- predvar[ind]
  #
  ################################################################################
  # SIMULATIONS
  #
  if(sim) {
    # SIMULATE COEFFICIENTS
    k <- length(coef)
    eigen <- eigen(vcov)
    X <- matrix(rnorm(length(coef)*nsim),nsim)
    coefsim <- coef + eigen$vectors %*% diag(sqrt(eigen$values),k) %*% t(X)
    # COMPUTE MINIMUM
    minsim <- apply(coefsim,2,function(coefi) {
      pred <- drop(Xpredall%*%coefi)
      ind <- which.max(pred)
      return(predvar[ind])
    })
  }
  #
  ################################################################################
  #
  res <- if(sim) minsim else min
  #
  return(res)
}



##-------------------------------------------------------------------------------------####
## 3. FUNCTION OF CASE-CROSSOVER ANALYSIS   lag-06 or other                            ####
##-------------------------------------------------------------------------------------####
#data  <- subset(cc2, per_mort == "Yes" & country == "Benin")
#lag <- 6
##confounder = con
#varper=c(10,75,90)
#lagnk=2
#cen=list(min=TRUE, pct=c(25,50,75))
#estpct=c(1,2.5,5,10,90,95,97.5,99)
#fml <- as.formula(paste0("status~cb_temp+strata(alert_id)+factor(",paste0(confounder,collapse = "+"), ")"))

casecrs <- function (status, id, confounder=NULL, lag, varper, lagnk, cen=list(min=NULL, max=NULL, degree=NULL, pct=NULL), estpct, data){
  
  ## Input: status, id, confounder: variables used to define the formula applied to clogit() in the form:
  ##                                case.status~exposure+confounder+strata(matched.set)
  ##                                status: case status, 1=case, 0=control
  ##                                id: ID for participants
  ##                                confounder: optional, vector of covariates to be included in the model
  ##        lag: the maximum lag in the cross basis
  ##        varper: numeric vector of percentiles of the distribution of temperature for internal knots
  ##        lagnk: the number of internal knots in the lag-response dimension
  ##        cen: a list to define the centering temperature
  ##            - "min" and "max": optional, "TRUE" if the minimum or maximum mortality temperature to be used
  ##            - "degree": optional, numeric vector of temperature (?C)
  ##            - "pct": optional, numeric vector of the percentiles of temperature distribution
  ##        estpct: numeric vector of percentiles of temperature distributions for effect estimate compared to centering temperature
  ## Build cross-basis function of temperature and lags
  ## Note: (1) exposure-response: natural cubic spline with internal knots placed at percentile of the temperature 
  ##                              distribution as defined by "varper"
  ##       (2) lag-response: natural cubic spline with an intercept and n="lag" internal knots placed at 
  ##                         equally spaced values on the log scale
  
  ## 1. delete observations with NA in temperature
  dat <- subset(data, rowSums(is.na(data[which(names(data)%in%paste0("temp_s",0:lag))]))==0)
  ## 2. extract matrix of temperature at lag0 to lag="lag"
  mat <- c("temp_s0", "temp_s1", "temp_s2", "temp_s3", "temp_s4", "temp_s5", "temp_s6")
  mat_temp <- as.matrix(dat[mat])
  
  # mat_temp <- as.matrix(select(dat, all_of(paste0("temp_s",0:lag))))
  
  ## 3. define basis for temperature
  argvar <- list(fun="ns", knots=quantile(mat_temp, varper/100, na.rm=T))
  ## 4. define basis for lag
  arglag <- list(fun="ns", knots=logknots(lag,lagnk))
  ## 5. build the cross-basis function
  cb_temp <- crossbasis(mat_temp, lag=c(0,lag), argvar=argvar, arglag=arglag)
  
  ## Different percentile of the temperature matrix
  tper <- quantile(mat_temp,seq(0,100,1)/100)
  
  ## Temperature summary for case days
  tsum_case       <- summary(subset(dat, status==1)$temp_s0)
  tsum_case["SD"] <- sd(subset(dat, status==1)$temp_s0)
  
  ## Temperature summary for control days
  tsum_control       <- summary(subset(dat, status==0)$temp_s0)
  tsum_control["SD"] <- sd(subset(dat, status==0)$temp_s0)
  
  ## Conditional logistic regression ####
  
  if (is.null(confounder)==F){
    fml <- as.formula(paste0(status,"~cb_temp+strata(",id,")+",paste0(confounder,collapse = "+")))
  } else {
    fml <- as.formula(paste0(status,"~cb_temp+strata(",id,")"))
  }
  mod <- clogit(fml,data=dat)
  # mod <- clogit(status~cb_temp+strata(alert_id),data=dat)
  # Reduction to overall cumulative (it is irrelevant the cen value)
  red <- crossreduce(cb_temp, mod, cen=20)
  # Store reduced coefs
  coef <- coef(red)
  vcov <- vcov(red)
  
  ## centering temperature
  cen_temp <- NULL;cen_name <- NULL
  
  if (is.null(cen$min)==F){
    cen_temp <- c(cen_temp, findmin(cb_temp, mod))
    cen_name <- c(cen_name,"min")
  }
  if (is.null(cen$max)==F){
    cen_temp <- c(cen_temp, findmax(cb_temp, mod))
    cen_name <- c(cen_name,"max")
  }
  if (is.null(cen$degree)==F){
    cen_temp <- c(cen_temp, cen$degree)
    cen_name <- c(cen_name, paste0(cen$degree," degree"))
  }
  if (is.null(cen$pct)==F){
    cen_temp <- c(cen_temp, quantile(mat_temp,cen$pct/100))
    cen_name <- c(cen_name, paste0(cen$pct,"th"))
  }
  
  ## Predict ORs from each cen_temp to each estpct
  estimate             <- list()
  for (i in 1:length(cen_temp)){
    pred               <- crosspred(cb_temp, mod, model.link="logit", cen=cen_temp[i], at=quantile(mat_temp,estpct/100))
    estimate[[i]]      <-  round(data.frame(beta=pred$allfit, SE=pred$allse), 3)
    estimate[[i]]$low       <-  estimate[[i]]$beta-qnorm(1-0.05/2)*estimate[[i]]$SE
    estimate[[i]]$high      <-  estimate[[i]]$beta+qnorm(1-0.05/2)*estimate[[i]]$SE
    estimate[[i]]$OR        <-  exp(estimate[[i]]$beta)
    estimate[[i]]$CIlow     <-  exp(estimate[[i]]$low)
    estimate[[i]]$CIhigh    <-  exp(estimate[[i]]$high)
    estimate[[i]]$temp <- as.numeric(rownames(estimate[[i]]))
    estimate[[i]]$perc <- paste0(estpct,"th")
    estimate[[i]]$cen  <- cen_name[i]
    estimate[[i]]$cen_min <- findmin(cb_temp, mod)
    estimate[[i]]$cen_50 <- quantile(mat_temp, 0.5, na.rm=TRUE)
    estimate[[i]]$cen_75 <- quantile(mat_temp, 0.75, na.rm=TRUE)
    estimate[[i]]      <- dplyr::select(estimate[[i]],c(cen,perc,temp,everything()))
  }
  estimate_all <- do.call(rbind,estimate)
  rownames(estimate_all) <- NULL
  
  ## output:result, a list containing the following elements
  ##        - n_case: number of cases
  ##        - n_control: number of controls
  ##        - tper: temperature distribution (percentiles) 
  ##        - tsum_case: summary of temperature on case days
  ##        - tsum_control: summary of temperature on case days
  ##        - coef: coefficients for the overall association
  ##        - vcov: variance-covariance of coefs for overall association
  ##        - estimate: OR and CI at the "estpct" percentile of temperature distribution compared to each cen_temp
  ##        output for plots
  ##        - mat_temp: matrix of temperature 
  ##        - cb_temp: cross-basis of temperature
  ##        - model_coef: coefficients of conditional logistic regression model
  ##        - model_vcov: variance matrix of conditional logistic regression model
  
  result              <- NULL
  result$n_case       <- nrow(subset(dat,status==1))
  result$n_control    <- nrow(subset(dat,status==0))
  result$tper         <- tper
  result$tsum_case    <- tsum_case
  result$tsum_control <- tsum_control
  result$coef         <- coef
  result$vcov         <- vcov
  result$estimate     <- estimate_all
  result$mat_temp     <- mat_temp
  result$cb_temp      <- cb_temp
  result$model_coef   <- mod$coefficients
  result$model_vcov   <- mod$var
  return(result)
}


#########################################################################################################
##### Stillbirths
#########################################################################################################

### Define list of provinces to be entered in the following loop
country=unique(alert$country)

### Generate empty lists to fill with first stage estimates
coef_list <- vector(mode = "list", length = length(country))
vcov_list <- vector(mode = "list", length = length(country))

### Define set of percentiles where I'll get the predictions
predper <- c(1:99)

### Generate an empty matrix to fill with province-specific temperature percentiles
tmeanper_matrix=matrix(data=NA, nrow=length(country), ncol=length(predper))

# start with firt vector
dlist <- list(cc.still.ben, cc.still.mal, cc.still.tan, cc.still.ug)

for (i in 1:length(dlist)){
  # i = 1
  data <- dlist[[i]]
  
  lag=6
  lagnk=2
  varper=c(10,75,90)
  cen=list(min=TRUE, pct=c(25,50,75))
  estpct=c(1,2.5,5,10,90,95,97.5,99)
  
  ## 1. delete observations with NA in temperature
  dat <- subset(data, rowSums(is.na(data[which(names(data)%in%paste0("temp_s",0:lag))]))==0)
  ## 2. extract matrix of temperature at lag0 to lag="lag"
  mat <- c("temp_s0", "temp_s1", "temp_s2", "temp_s3", "temp_s4", "temp_s5", "temp_s6")
  mat_temp <- as.matrix(dat[mat])
  
  # mat_temp <- as.matrix(select(dat, all_of(paste0("temp_s",0:lag))))
  
  ## 3. define basis for temperature
  argvar <- list(fun="ns", knots=quantile(mat_temp, varper/100, na.rm=T))
  ## 4. define basis for lag
  arglag <- list(fun="ns", knots=logknots(lag,lagnk))
  ## 5. build the cross-basis function
  cb_temp <- crossbasis(mat_temp, lag=c(0,lag), argvar=argvar, arglag=arglag)
  
  ## Different percentile of the temperature matrix
  tper <- quantile(mat_temp,seq(1,99,1)/100)
  
  ## Temperature summary for case days
  tsum_case       <- summary(subset(dat, status==1)$temp_s0)
  tsum_case["SD"] <- sd(subset(dat, status==1)$temp_s0)
  
  ## Temperature summary for control days
  tsum_control       <- summary(subset(dat, status==0)$temp_s0)
  tsum_control["SD"] <- sd(subset(dat, status==0)$temp_s0)
  
  ## Conditional logistic regression 
  mod <- clogit(status~cb_temp+strata(alert_id),data=dat)
  # Reduction to overall cumulative (it is irrelevant the cen value)
  red <- crossreduce(cb_temp, mod, cen=20)
  
  # Store reduced coefs
  coef_list[[i]]=coef(red)
  vcov_list[[i]]=vcov(red)
  
  tmean_prov=quantile(dat$temp_s0, na.rm=T, predper/100)
  tmeanper_matrix[i,]=tmean_prov
  
  print(i)
}

######### meta-regression
### Extract the parameters
coef_df <- t(sapply(coef_list, function(x) x))
vcov_df <- lapply(vcov_list, function(x) x)

### fit a multivariate meta-analysis model
metareg <- mixmeta(coef_df ~ 1, vcov_df, method="ml", control=list(showiter=F))
summary(metareg)

### extract distribution of temperatures
tmeanper <- colMeans(tmeanper_matrix)
names(tmeanper)=names(tmean_prov)
per <- names(tmeanper)

### define the exposure response function
bvar <- crossbasis(tmeanper, argvar=argvar, arglag=arglag, Boundary.knots=tmeanper[c("1%","99%")])

### predict pooled estimates at 75th versus 99 (main approach)
cpmeta.75_still  <- crosspred(bvar, coef=coef(metareg), vcov=vcov(metareg), model.link="log", at=tmeanper, cen=tmeanper[c("75%")])
cpmeta.90_still  <- crosspred(bvar, coef=coef(metareg), vcov=vcov(metareg), model.link="log", at=tmeanper, cen=tmeanper[c("90%")])

### predict pooled estimates at MMT versus 99 (sensitivity approach)
cen        <- tmeanper["50%"]
cpmeta     <- crosspred(bvar, coef=coef(metareg), vcov=vcov(metareg), model.link="log", at=tmeanper, cen=cen)
mmt=cpmeta$predvar[which(cpmeta$allRRfit==min(cpmeta$allRRfit))]

# SPECIFIC FOR EACH OUTCOME - STILLBIRTH
cpmeta.MMT_still <- crosspred(bvar, coef=coef(metareg), vcov=vcov(metareg), model.link="log", at=tmeanper, cen=mmt)
pp <- names(tmeanper[tmeanper == mmt])
pp_still <- as.numeric(gsub('%','',pp))
tmeanper_still <- tmeanper

#########################################################################################################
##### INTRAPARTUM STILLBIRTH
#########################################################################################################

### 1. Effect estimation (refer to the function "casecrs" in "CHAMNHA_functions.R") 
### Define list of provinces to be entered in the following loop
country=unique(alert$country)

### Generate empty lists to fill with first stage estimates
coef_list <- vector(mode = "list", length = length(country))
vcov_list <- vector(mode = "list", length = length(country))

### Define set of percentiles where I'll get the predictions
#predper <- c(seq(0,1,0.1), 2:98, seq(99,100,0.1))
predper <- c(1:99)
### Generate an empty matrix to fill with province-specific temperature percentiles
tmeanper_matrix=matrix(data=NA, nrow=length(country), ncol=length(predper))

# start with firt vector
dlist <- list( cc.intra.ben, cc.intra.mal, cc.intra.tan, cc.intra.ug)


for (i in 1:length(dlist)){
  # i = 1
  data <- dlist[[i]]
  
  lag=6
  lagnk=2
  varper=c(10,75,90)
  cen=list(min=TRUE, pct=c(25,50,75))
  estpct=c(1,2.5,5,10,90,95,97.5,99)
  
  ## 1. delete observations with NA in temperature
  dat <- subset(data, rowSums(is.na(data[which(names(data)%in%paste0("temp_s",0:lag))]))==0)
  ## 2. extract matrix of temperature at lag0 to lag="lag"
  mat <- c("temp_s0", "temp_s1", "temp_s2", "temp_s3", "temp_s4", "temp_s5", "temp_s6")
  mat_temp <- as.matrix(dat[mat])
  
  # mat_temp <- as.matrix(select(dat, all_of(paste0("temp_s",0:lag))))
  
  ## 3. define basis for temperature
  argvar <- list(fun="ns", knots=quantile(mat_temp, varper/100, na.rm=T))
  ## 4. define basis for lag
  arglag <- list(fun="ns", knots=logknots(lag,lagnk))
  ## 5. build the cross-basis function
  cb_temp <- crossbasis(mat_temp, lag=c(0,lag), argvar=argvar, arglag=arglag)
  
  ## Different percentile of the temperature matrix
  tper <- quantile(mat_temp,seq(1,99,1)/100)
  
  ## Temperature summary for case days
  tsum_case       <- summary(subset(dat, status==1)$temp_s0)
  tsum_case["SD"] <- sd(subset(dat, status==1)$temp_s0)
  
  ## Temperature summary for control days
  tsum_control       <- summary(subset(dat, status==0)$temp_s0)
  tsum_control["SD"] <- sd(subset(dat, status==0)$temp_s0)
  
  ## Conditional logistic regression 
  mod <- clogit(status~cb_temp+strata(alert_id),data=dat)
  # Reduction to overall cumulative (it is irrelevant the cen value)
  red <- crossreduce(cb_temp, mod, cen=20)
  
  # Store reduced coefs
  coef <- coef(red)
  vcov <- vcov(red)
  
  coef_list[[i]]=coef(red)
  vcov_list[[i]]=vcov(red)
  
  tmean_prov=quantile(dat$temp_s0, na.rm=T, predper/100)
  tmeanper_matrix[i,]=tmean_prov
  
  print(i)
}

######### meta-regression
### Extract the parameters
coef_df <- t(sapply(coef_list, function(x) x))
vcov_df <- lapply(vcov_list, function(x) x)

### fit a multivariate meta-analysis model
metareg <- mixmeta(coef_df ~ 1, vcov_df, method="ml", control=list(showiter=F))
summary(metareg)

### extract distribution of temperatures
tmeanper <- colMeans(tmeanper_matrix)
names(tmeanper)=names(tmean_prov)
per <- names(tmeanper)

### define the exposure response function
bvar <- crossbasis(tmeanper, argvar=argvar, arglag=arglag, Boundary.knots=tmeanper[c("1%","99%")])

### Compute the BLUPS
# blups <- blup(metareg,vcov=T)

### predict pooled estimates at 75th versus 99 (main approach)
cpmeta.75_fresh  <- crosspred(bvar, coef=coef(metareg), vcov=vcov(metareg), model.link="log", at=tmeanper, cen=tmeanper[c("75%")])
cpmeta.90_fresh  <- crosspred(bvar, coef=coef(metareg), vcov=vcov(metareg), model.link="log", at=tmeanper, cen=tmeanper[c("90%")])

### predict pooled estimates at MMT versus 99 (sensitivity approach)
cen        <- tmeanper["50%"]
cpmeta     <- crosspred(bvar, coef=coef(metareg), vcov=vcov(metareg), model.link="log", at=tmeanper, cen=cen)
mmt=cpmeta$predvar[which(cpmeta$allRRfit==min(cpmeta$allRRfit))]

# SPECIFIC FOR EACH OUTCOME - Intrapartum stillbirths
cpmeta.MMT_fresh <- crosspred(bvar, coef=coef(metareg), vcov=vcov(metareg), model.link="log", at=tmeanper, cen=mmt)
pp <- names(tmeanper[tmeanper == mmt])
pp_fresh <- as.numeric(gsub('%','',pp))
tmeanper_fresh <- tmeanper


#########################################################################################################
##### ANTERPARTUM STILLBIRTH
#########################################################################################################

### Define list of provinces to be entered in the following loop
country=unique(alert$country)

### Generate empty lists to fill with first stage estimates
coef_list <- vector(mode = "list", length = length(country))
vcov_list <- vector(mode = "list", length = length(country))

### Define set of percentiles where I'll get the predictions
predper <- c(1:99)

### Generate an empty matrix to fill with province-specific temperature percentiles
tmeanper_matrix=matrix(data=NA, nrow=length(country), ncol=length(predper))

# start with firt vector
dlist <- list( cc.ante.ben, cc.ante.mal, cc.ante.tan, cc.ante.ug)


for (i in 1:length(dlist)){
  # i = 1
  data <- dlist[[i]]
  
  lag=6
  lagnk=2
  varper=c(10,75,90)
  cen=list(min=TRUE, pct=c(25,50,75))
  estpct=c(1,2.5,5,10,90,95,97.5,99)
  
  ## 1. delete observations with NA in temperature
  dat <- subset(data, rowSums(is.na(data[which(names(data)%in%paste0("temp_s",0:lag))]))==0)
  ## 2. extract matrix of temperature at lag0 to lag="lag"
  mat <- c("temp_s0", "temp_s1", "temp_s2", "temp_s3", "temp_s4", "temp_s5", "temp_s6")
  mat_temp <- as.matrix(dat[mat])
  
  # mat_temp <- as.matrix(select(dat, all_of(paste0("temp_s",0:lag))))
  
  ## 3. define basis for temperature
  argvar <- list(fun="ns", knots=quantile(mat_temp, varper/100, na.rm=T))
  ## 4. define basis for lag
  arglag <- list(fun="ns", knots=logknots(lag,lagnk))
  ## 5. build the cross-basis function
  cb_temp <- crossbasis(mat_temp, lag=c(0,lag), argvar=argvar, arglag=arglag)
  
  ## Different percentile of the temperature matrix
  tper <- quantile(mat_temp,seq(1,99,1)/100)
  
  ## Temperature summary for case days
  tsum_case       <- summary(subset(dat, status==1)$temp_s0)
  tsum_case["SD"] <- sd(subset(dat, status==1)$temp_s0)
  
  ## Temperature summary for control days
  tsum_control       <- summary(subset(dat, status==0)$temp_s0)
  tsum_control["SD"] <- sd(subset(dat, status==0)$temp_s0)
  
  ## Conditional logistic regression 
  mod <- clogit(status~cb_temp+strata(alert_id),data=dat)
  # Reduction to overall cumulative (it is irrelevant the cen value)
  red <- crossreduce(cb_temp, mod, cen=20)
  
  # Store reduced coefs
  coef <- coef(red)
  vcov <- vcov(red)
  
  coef_list[[i]]=coef(red)
  vcov_list[[i]]=vcov(red)
  
  tmean_prov=quantile(dat$temp_s0, na.rm=T, predper/100)
  tmeanper_matrix[i,]=tmean_prov
  
  print(i)
}

######### meta-regression
### Extract the parameters
coef_df <- t(sapply(coef_list, function(x) x))
vcov_df <- lapply(vcov_list, function(x) x)

### fit a multivariate meta-analysis model
metareg <- mixmeta(coef_df ~ 1, vcov_df, method="ml", control=list(showiter=F))
summary(metareg)

### extract distribution of temperatures
tmeanper <- colMeans(tmeanper_matrix)
names(tmeanper)=names(tmean_prov)
per <- names(tmeanper)

### define the exposure response function
bvar <- crossbasis(tmeanper, argvar=argvar, arglag=arglag, Boundary.knots=tmeanper[c("1%","99%")])

### predict pooled estimates at 75th versus 99 (main approach)
cpmeta.75_mac  <- crosspred(bvar, coef=coef(metareg), vcov=vcov(metareg), model.link="log", at=tmeanper, cen=tmeanper[c("75%")])
cpmeta.90_mac  <- crosspred(bvar, coef=coef(metareg), vcov=vcov(metareg), model.link="log", at=tmeanper, cen=tmeanper[c("90%")])

### predict pooled estimates at MMT versus 99 (sensitivity approach)
cen        <- tmeanper["50%"]
cpmeta     <- crosspred(bvar, coef=coef(metareg), vcov=vcov(metareg), model.link="log", at=tmeanper, cen=cen)
mmt=cpmeta$predvar[which(cpmeta$allRRfit==min(cpmeta$allRRfit))]

# SPECIFIC FOR EACH OUTCOME - Antepartum STILLBIRTHS
cpmeta.MMT_mac <- crosspred(bvar, coef=coef(metareg), vcov=vcov(metareg), model.link="log", at=tmeanper, cen=mmt)
pp <- names(tmeanper[tmeanper == mmt])
pp_mac <- as.numeric(gsub('%','',pp))
tmeanper_mac <- tmeanper

#########################################################################################################
##### PERINATAL DEATHS
#########################################################################################################

### Define list of provinces to be entered in the following loop
country=unique(alert$country)

### Generate empty lists to fill with first stage estimates
coef_list <- vector(mode = "list", length = length(country))
vcov_list <- vector(mode = "list", length = length(country))

### Define set of percentiles where I'll get the predictions
#predper <- c(seq(0,1,0.1), 2:98, seq(99,100,0.1))
predper <- c(1:99)
### Generate an empty matrix to fill with province-specific temperature percentiles
tmeanper_matrix=matrix(data=NA, nrow=length(country), ncol=length(predper))

# start with firt vector
dlist <- list( cc.per.ben, cc.per.mal, cc.per.tan, cc.per.ug)


for (i in 1:length(dlist)){
  # i = 1
  data <- dlist[[i]]
  
  lag=6
  lagnk=2
  varper=c(10,75,90)
  cen=list(min=TRUE, pct=c(25,50,75))
  estpct=c(1,2.5,5,10,90,95,97.5,99)
  
  ## 1. delete observations with NA in temperature
  dat <- subset(data, rowSums(is.na(data[which(names(data)%in%paste0("temp_s",0:lag))]))==0)
  ## 2. extract matrix of temperature at lag0 to lag="lag"
  mat <- c("temp_s0", "temp_s1", "temp_s2", "temp_s3", "temp_s4", "temp_s5", "temp_s6")
  mat_temp <- as.matrix(dat[mat])
  
  # mat_temp <- as.matrix(select(dat, all_of(paste0("temp_s",0:lag))))
  
  ## 3. define basis for temperature
  argvar <- list(fun="ns", knots=quantile(mat_temp, varper/100, na.rm=T))
  ## 4. define basis for lag
  arglag <- list(fun="ns", knots=logknots(lag,lagnk))
  ## 5. build the cross-basis function
  cb_temp <- crossbasis(mat_temp, lag=c(0,lag), argvar=argvar, arglag=arglag)
  
  ## Different percentile of the temperature matrix
  tper <- quantile(mat_temp,seq(1,99,1)/100)
  
  ## Temperature summary for case days
  tsum_case       <- summary(subset(dat, status==1)$temp_s0)
  tsum_case["SD"] <- sd(subset(dat, status==1)$temp_s0)
  
  ## Temperature summary for control days
  tsum_control       <- summary(subset(dat, status==0)$temp_s0)
  tsum_control["SD"] <- sd(subset(dat, status==0)$temp_s0)
  
  ## Conditional logistic regression 
  mod <- clogit(status~cb_temp+strata(alert_id),data=dat)
  # Reduction to overall cumulative (it is irrelevant the cen value)
  red <- crossreduce(cb_temp, mod, cen=20)
  
  # Store reduced coefs
  coef <- coef(red)
  vcov <- vcov(red)
  
  coef_list[[i]]=coef(red)
  vcov_list[[i]]=vcov(red)
  
  tmean_prov=quantile(dat$temp_s0, na.rm=T, predper/100)
  tmeanper_matrix[i,]=tmean_prov
  
  print(i)
}

######### meta-regression
### Extract the parameters
coef_df <- t(sapply(coef_list, function(x) x))
vcov_df <- lapply(vcov_list, function(x) x)

### fit a multivariate meta-analysis model
metareg <- mixmeta(coef_df ~ 1, vcov_df, method="ml", control=list(showiter=F))
summary(metareg)

### extract distribution of temperatures
tmeanper <- colMeans(tmeanper_matrix)
names(tmeanper)=names(tmean_prov)
per <- names(tmeanper)

### define the exposure response function
bvar <- crossbasis(tmeanper, argvar=argvar, arglag=arglag, Boundary.knots=tmeanper[c("1%","99%")])

### predict pooled estimates at 75th versus 99 (main approach)
cpmeta.75_per  <- crosspred(bvar, coef=coef(metareg), vcov=vcov(metareg), model.link="log", at=tmeanper, cen=tmeanper[c("75%")])
cpmeta.90_per  <- crosspred(bvar, coef=coef(metareg), vcov=vcov(metareg), model.link="log", at=tmeanper, cen=tmeanper[c("90%")])

### predict pooled estimates at MMT versus 99 (sensitivity approach)
cen        <- tmeanper["50%"]
cpmeta     <- crosspred(bvar, coef=coef(metareg), vcov=vcov(metareg), model.link="log", at=tmeanper, cen=cen)
mmt=cpmeta$predvar[which(cpmeta$allRRfit==min(cpmeta$allRRfit))]

# SPECIFIC FOR EACH OUTCOME - perinatal deaths
cpmeta.MMT_per <- crosspred(bvar, coef=coef(metareg), vcov=vcov(metareg), model.link="log", at=tmeanper, cen=mmt)
pp <- names(tmeanper[tmeanper == mmt])
pp_per <- as.numeric(gsub('%','',pp))
tmeanper_per <- tmeanper


#########################################################################################################
##### CREATE FIGURE 3 of manuscript
#########################################################################################################

png("05_figures/mean_temp_meta-regression_cen75_9x9.png", width = 9, height = 9, units = "in", res = 300)
par(mfrow=c(2,2), mar=c(4, 4, 2, 1) + 0.1, oma=c(0, 0, 2, 0))

# stillbirths
plot(cpmeta.75_still, "overall", lwd=1.5,  ylim=c(0.5, 5), col=2, axes=F, xlab="", ylab="")
mtext("A) Stillbirths", cex=1.0, line=1.0, font=1, family = "sans")  
axis(side=1, at=seq(0, 40,1.0), cex.axis=0.9)
mtext(side=1, "Air temperature (°C)", font=1, line=2.5, cex=0.8, family = "sans")
axis(side=2,  at=seq(0.5,5.0,0.5), cex.axis=0.9)
mtext(side=2, "RR (95% CI)", font=1, line=2.5, cex=0.8, family = "sans")


# Antepartum
plot(cpmeta.75_mac, "overall", lwd=1.5,  ylim=c(0.5, 5), col=2, axes=F, xlab="", ylab="")
mtext("B) Antepartum stillbirths", cex=1.0, line=1.0, font=1, family = "sans")   
axis(side=1, at=seq(0, 40,1.0), cex.axis=0.9)
mtext(side=1, "Air temperature (°C)", font=1, line=2.5, cex=0.8, family = "sans")
axis(side=2,  at=seq(0.5,5.0,0.5), cex.axis=0.9)
mtext(side=2, "RR (95% CI)", font=1, line=2.5, cex=0.8, family = "sans")

# Intrapartum stillbirths
plot(cpmeta.75_fresh, "overall", lwd=1.5,  ylim=c(0.5, 5), col=2, axes=F, xlab="", ylab="")
mtext("C) Intrapartum stillbirths", cex=1.0, line=1.0, font=1, family = "sans")  
axis(side=1, at=seq(0, 40,1.0), cex.axis=0.9)
mtext(side=1, "Air temperature (°C)", font=1, line=2.5, cex=0.8, family = "sans")
axis(side=2,  at=seq(0.5,5.0,0.5), cex.axis=0.9)
mtext(side=2, "RR (95% CI)", font=1, line=2.5, cex=0.8, family = "sans")

# PERINTATAL
plot(cpmeta.75_per, "overall", lwd=1.5,  ylim=c(0.5, 5), col=2, axes=F, xlab="", ylab="")
mtext("D) Perinatal deaths", cex=1.0, line=1.0, font=1, family = "sans") 
axis(side=1, at=seq(0, 40,1.0), cex.axis=0.9)
mtext(side=1, "Air temperature (°C)", font=1, line=2.5, cex=0.8, family = "sans")
axis(side=2,  at=seq(0.5,5.0,0.5), cex.axis=0.9)
mtext(side=2, "RR (95% CI)", font=1, line=2.5, cex=0.8, family = "sans")

dev.off()

