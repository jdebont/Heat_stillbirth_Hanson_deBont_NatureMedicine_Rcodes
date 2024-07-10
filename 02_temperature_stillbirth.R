##############################################################################
#
# Contact: Jeroen de Bont (jeroen.de.bont@ki.se)
#
##############################################################################
#
# Script:  Analyse the association between mean temperature (lag-06) and stillbirths
#
##############################################################################

#add the necessary libraries
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


#####################################################
##### analyses by country and stillbirth outcome
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

## create a list containing all datasets
dlist <- list(cc.still.ben, cc.intra.ben, cc.ante.ben, cc.per.ben,
              cc.still.mal, cc.intra.mal, cc.ante.mal, cc.per.mal,
              cc.still.tan, cc.intra.tan, cc.ante.tan, cc.per.tan,
              cc.still.ug, cc.intra.ug, cc.ante.ug, cc.per.ug)

## names of the outcomes in the order of that in the dlist
name_outcome <- c("still.ben","intra.ben","ante.ben", "per.ben",
                  "still.mal","intra.mal","ante.mal", "per.mal",
                  "still.tan","intra.tan","ante.tan", "per.tan",
                  "still.ug","intra.ug","ante.ug", "per.ug")


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
    cdf <- ecdf(mat_temp)
    estimate[[i]]$cen_min_perc <- cdf(estimate[[i]]$cen_min)*100
    #estimate[[i]]$cen_min_perc_doub <- quantile(mat_temp, cen_min_perc/100, na.rm=TRUE)
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

##-------------------------------------------------------------------------------------####
## 4. FUNCTION OF CASE-CROSSOVER ANALYSIS   lag-0                                      ####
##-------------------------------------------------------------------------------------####

casecrs_0 <- function (status, id, confounder=NULL, lag, varper, lagnk, cen=list(min=NULL, max=NULL, degree=NULL, pct=NULL), estpct, data){
  
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
  #mat <- c("temp_s0", "temp_s1")
  #mat_temp <- as.matrix(dat[mat])
  mat_temp <- as.matrix(select(dat, all_of(paste0("temp_s",0:lag))))
  
  ## 3. define basis for temperature
  argvar <- list(fun="ns", knots=quantile(mat_temp, varper/100, na.rm=T))
  
  ## 4. define basis for lag
  arglag <- list(fun="ns")
  
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
  
  # if (is.null(confounder)==F){
  #   fml <- as.formula(paste0(status,"~cb_temp+strata(",id,")+",paste0(confounder,collapse = "+")))
  # } else {
  #   fml <- as.formula(paste0(status,"~cb_temp+strata(",id,")"))
  # }
  mod <- clogit(status~cb_temp+strata(alert_id), data=dat)
  
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
##### Run the models
######################################################################################################### 

##Run the analyses after creating the specified functions above
### 1. Effect estimation ####
result1 <- list()
for (i in 1:length(dlist)){
  # i = 1
  dat <- dlist[[i]]
  dat <- as.data.frame(dat)
  result1[[i]] <- casecrs(status="status", id="alert_id",  lag=6, varper=c(10,75,90), lagnk=2, cen=list(min=TRUE, pct=c(25,50,75,90)), estpct=c(1,2.5,5,10,90,95,97.5,99), data=dat)
  print(i)
}

names(result1) <- name_outcome

#### Put together the main results
still.ben  <- cbind(as.data.table(result1[["still.ben"]]$estimate), model = "All stillbirth", country = "Benin")
intra.ben  <- cbind(as.data.table(result1[["intra.ben"]]$estimate), model = "Intra stillbirth", country = "Benin")
ante.ben   <- cbind(as.data.table(result1[["ante.ben"]]$estimate), model = "Antermortum", country = "Benin")
per.ben    <- cbind(as.data.table(result1[["per.ben"]]$estimate), model = "Pernatal death", country = "Benin")

still.mal  <- cbind(as.data.table(result1[["still.mal"]]$estimate), model = "All stillbirth", country = "Malawi")
intra.mal  <- cbind(as.data.table(result1[["intra.mal"]]$estimate), model = "Intra stillbirth", country = "Malawi")
ante.mal   <- cbind(as.data.table(result1[["ante.mal"]]$estimate), model = "Antermortum", country = "Malawi")
per.mal    <- cbind(as.data.table(result1[["per.mal"]]$estimate), model = "Pernatal death", country = "Malawi")

still.tan  <- cbind(as.data.table(result1[["still.tan"]]$estimate), model = "All stillbirth", country = "Tanzania")
intra.tan  <- cbind(as.data.table(result1[["intra.tan"]]$estimate), model = "Intra stillbirth", country = "Tanzania")
ante.tan   <- cbind(as.data.table(result1[["ante.tan"]]$estimate), model = "Antermortum", country = "Tanzania")
per.tan    <- cbind(as.data.table(result1[["per.tan"]]$estimate), model = "Pernatal death", country = "Tanzania")

still.ug  <- cbind(as.data.table(result1[["still.ug"]]$estimate), model = "All stillbirth", country = "Uganda")
intra.ug  <- cbind(as.data.table(result1[["intra.ug"]]$estimate), model = "Intra stillbirth", country = "Uganda")
ante.ug   <- cbind(as.data.table(result1[["ante.ug"]]$estimate), model = "Antermortum", country = "Uganda")
per.ug    <- cbind(as.data.table(result1[["per.ug"]]$estimate), model = "Pernatal death", country = "Uganda")

all <- rbind(still.ben,intra.ben,ante.ben, per.ben,
             still.mal,intra.mal,ante.mal, per.mal,
             still.tan,intra.tan,ante.tan, per.tan,
             still.ug,intra.ug,ante.ug, per.ug)

# store all results
save(all, file = "results_country_still_9x9.Rda")


#########################################################################################################
##### Meta-analyse the results
######################################################################################################### 

# load all previous results from the other estiamtes and combine them
load("results_country_still_9x9.Rda")

# meta analyses
library(meta)

###### CREATE META-ANALYSES RESULTS
sum <- all[ ( perc == "99th")] # use 99th percentile as our increase

## Create table 4: descriptive PCA analyses
model_st <- c("All stillbirth", "Intra stillbirth", "Antermortum", "Pernatal death")
cen <- c("min","25th", "50th", "75th") # we have estimated different "cen" values. From "min", 25th, 50th and 75th

# Pooled Estimate using a function by outcome and by centred value. This is used to estimate figure 2
model <- data.frame()  
for (i in model_st){ 
  for (j in cen){  
    
    # Pooled Estimate
    cohort <- sum[model == paste0(i) & cen == paste0(j),.(cen, perc, OR, CIlow, CIhigh, model, country)]
    
    meta_mmt <- metagen(log(OR), lower = log(CIlow), upper = log(CIhigh),
                        studlab = country, data=cohort, sm = "OR")
    
    # forest(meta_mmt)
    # extract
    res <-data.table( beta=meta_mmt$TE.random, se=meta_mmt$seTE.random, i2=meta_mmt$I2)
    res[, model := i]
    res[, cen := j]
    model <- as.data.table(rbind(model, res))
  }
}

### make table
model$low<-model$beta-qnorm(1-0.05/2)*model$se
model$high<-model$beta+qnorm(1-0.05/2)*model$se
model$or<-exp(model$beta)
model$or_low<-exp(model$low)
model$or_high<-exp(model$high)

# reduce digits and combine
model[, or := formatC(or,digits=2, format="f")]
model[, or_low := formatC(or_low,digits=2, format="f")]
model[, or_high := formatC(or_high, digits=2, format="f")]
model[,c('hr_all'):=paste0(or," (", or_low, "; ",or_high,")")]

# add I2
model[, i2 := formatC(as.numeric(i2),digits=2, format="f")]
save(model, file = "04_results/tmean_lag06_9x9.Rda")


#########################################################################################################
##### Create figure 2
######################################################################################################### 

# Pooled Estimate

### Stillbirths

cohort <- sum[model == "All stillbirth" & cen == "75th",.(cen, perc, OR, CIlow, CIhigh, model, country)]
meta_stil <- metagen(log(OR), lower = log(CIlow), upper = log(CIhigh),
                     studlab = country, data=cohort, sm = "OR")

cohort <- sum[model == "Intra stillbirth" & cen == "75th",.(cen, perc, OR, CIlow, CIhigh, model, country)]
meta_fre <- metagen(log(OR), lower = log(CIlow), upper = log(CIhigh),
                    studlab = country, data=cohort, sm = "OR")

cohort <- sum[model == "Antermortum" & cen == "75th",.(cen, perc, OR, CIlow, CIhigh, model, country)]
meta_mac <- metagen(log(OR), lower = log(CIlow), upper = log(CIhigh),
                    studlab = country, data=cohort, sm = "OR")

cohort <- sum[model == "Pernatal death" & cen == "75th",.(cen, perc, OR, CIlow, CIhigh, model, country)]
meta_per <- metagen(log(OR), lower = log(CIlow), upper = log(CIhigh),
                    studlab = country, data=cohort, sm = "OR")


tiff("05_figures/meta_still_9x9.tiff", width = 1750, height = 650, res = 300)
forest(meta_stil, fixed = F, print.tau2 = F, print.pval.Q = F, 
       leftcols = c("country"), leftlabs = c("Country"), 
       rightcols = c("effect.ci"))

dev.off()

tiff("05_figures/meta_mac_9x9.tiff", width = 1750, height = 650, res = 300)
forest(meta_mac, fixed = F, print.tau2 = F, print.pval.Q = F,
       leftcols = c("country"), leftlabs = c("Country"),
       rightcols = c("effect.ci"))
dev.off()

tiff("05_figures/meta_fre_9x9.tiff", width = 1750, height = 650, res = 300)
forest(meta_fre, fixed = F, print.tau2 = F, print.pval.Q = F,
       leftcols = c("country"), leftlabs = c("Country"),
       rightcols = c("effect.ci"))
dev.off()

tiff("05_figures/meta_per_9x9.tiff", width = 1750, height = 650, res = 300)
forest(meta_per, fixed = F, print.tau2 = F, print.pval.Q = F,
       leftcols = c("country"), leftlabs = c("Country"),
       rightcols = c("effect.ci"))
dev.off()



#########################################################################################################
##### Associations by other lags
######################################################################################################### 

#########################################################################################################
##### Lag 0
######################################################################################################### 

### 1. Effect estimation using the casecrs__0 specific for lag 0
result1 <- list()
for (i in 1:length(dlist)){
  # i = 1
  dat <- dlist[[i]]
  dat <- as.data.frame(dat)
  result1[[i]] <- casecrs_0(status="status", id="alert_id",  lag=0, varper=c(10,75,90), lagnk=2, cen=list(min=TRUE, pct=c(25,50,75)), estpct=c(1,2.5,5,10,90,95,97.5,99), data=dat)
  print(i)
}

names(result1) <- name_outcome

#### Put together the main results
still.ben  <- cbind(as.data.table(result1[["still.ben"]]$estimate), model = "All stillbirth", country = "Benin")
intra.ben  <- cbind(as.data.table(result1[["intra.ben"]]$estimate), model = "Intra stillbirth", country = "Benin")
ante.ben   <- cbind(as.data.table(result1[["ante.ben"]]$estimate), model = "Antermortum", country = "Benin")
per.ben    <- cbind(as.data.table(result1[["per.ben"]]$estimate), model = "Pernatal death", country = "Benin")

still.mal  <- cbind(as.data.table(result1[["still.mal"]]$estimate), model = "All stillbirth", country = "Malawi")
intra.mal  <- cbind(as.data.table(result1[["intra.mal"]]$estimate), model = "Intra stillbirth", country = "Malawi")
ante.mal   <- cbind(as.data.table(result1[["ante.mal"]]$estimate), model = "Antermortum", country = "Malawi")
per.mal    <- cbind(as.data.table(result1[["per.mal"]]$estimate), model = "Pernatal death", country = "Malawi")

still.tan  <- cbind(as.data.table(result1[["still.tan"]]$estimate), model = "All stillbirth", country = "Tanzania")
intra.tan  <- cbind(as.data.table(result1[["intra.tan"]]$estimate), model = "Intra stillbirth", country = "Tanzania")
ante.tan   <- cbind(as.data.table(result1[["ante.tan"]]$estimate), model = "Antermortum", country = "Tanzania")
per.tan    <- cbind(as.data.table(result1[["per.tan"]]$estimate), model = "Pernatal death", country = "Tanzania")

still.ug  <- cbind(as.data.table(result1[["still.ug"]]$estimate), model = "All stillbirth", country = "Uganda")
intra.ug  <- cbind(as.data.table(result1[["intra.ug"]]$estimate), model = "Intra stillbirth", country = "Uganda")
ante.ug   <- cbind(as.data.table(result1[["ante.ug"]]$estimate), model = "Antermortum", country = "Uganda")
per.ug    <- cbind(as.data.table(result1[["per.ug"]]$estimate), model = "Pernatal death", country = "Uganda")

all <- rbind(still.ben,intra.ben,ante.ben, per.ben,
             still.mal,intra.mal,ante.mal,per.mal,
             still.tan,intra.tan,ante.tan, per.tan,
             still.ug,intra.ug,ante.ug, per.ug)

###### CREATE META-ANALYSES RESULTS
sum <- all[ ( perc == "99th")]

## Create table 4: descriptive PCA analyses
model_st <- c("All stillbirth", "Intra stillbirth", "Antermortum", "Pernatal death")
cen <- c("min","25th", "50th", "75th")

# Pooled Estimate
model <- data.frame()  
for (i in model_st){ 
  for (j in cen){  
    
    #i == "All stillbirth"
    #j == "min"  
    
    # Pooled Estimate
    cohort <- sum[model == paste0(i) & cen == paste0(j),.(cen, perc, OR, CIlow, CIhigh, model, country)]
    
    meta_mmt <- metagen(log(OR), lower = log(CIlow), upper = log(CIhigh),
                        studlab = country, data=cohort, sm = "OR")
    
    # extract
    res <-data.table( beta=meta_mmt$TE.random, se=meta_mmt$seTE.random, i2=meta_mmt$I2)
    res[, model := i]
    res[, cen := j]
    model <- as.data.table(rbind(model, res))
  }
}

### make table
model$low<-model$beta-qnorm(1-0.05/2)*model$se
model$high<-model$beta+qnorm(1-0.05/2)*model$se
model$or<-exp(model$beta)
model$or_low<-exp(model$low)
model$or_high<-exp(model$high)


model[, or := formatC(or,digits=2, format="f")]
model[, or_low := formatC(or_low,digits=2, format="f")]
model[, or_high := formatC(or_high, digits=2, format="f")]
model[,c('hr_all'):=paste0(or," (", or_low, "; ",or_high,")")]

# add I2
model[, i2 := formatC(as.numeric(i2),digits=2, format="f")]
save(model, file = "04_results/tmean_lag0_9x9.Rda")



#########################################################################################################
##### Lag 0-1
######################################################################################################### 


### 1. Effect estimation (refer to the function "casecrs" in "CHAMNHA_functions.R") ####
result1 <- list()
for (i in 1:length(dlist)){
  # i = 1
  dat <- dlist[[i]]
  dat <- as.data.frame(dat)
  result1[[i]] <- casecrs_0(status="status", id="alert_id",  lag=1, varper=c(10,75,90), lagnk=2, cen=list(min=TRUE, pct=c(25,50,75)), estpct=c(1,2.5,5,10,90,95,97.5,99), data=dat)
  print(i)
}

names(result1) <- name_outcome

#### Put together the main results
still.ben  <- cbind(as.data.table(result1[["still.ben"]]$estimate), model = "All stillbirth", country = "Benin")
intra.ben  <- cbind(as.data.table(result1[["intra.ben"]]$estimate), model = "Intra stillbirth", country = "Benin")
ante.ben   <- cbind(as.data.table(result1[["ante.ben"]]$estimate), model = "Antermortum", country = "Benin")
neo.ben    <- cbind(as.data.table(result1[["neo.ben"]]$estimate), model = "Neonatal death", country = "Benin")
per.ben    <- cbind(as.data.table(result1[["per.ben"]]$estimate), model = "Pernatal death", country = "Benin")

still.mal  <- cbind(as.data.table(result1[["still.mal"]]$estimate), model = "All stillbirth", country = "Malawi")
intra.mal  <- cbind(as.data.table(result1[["intra.mal"]]$estimate), model = "Intra stillbirth", country = "Malawi")
ante.mal   <- cbind(as.data.table(result1[["ante.mal"]]$estimate), model = "Antermortum", country = "Malawi")
neo.mal    <- cbind(as.data.table(result1[["neo.mal"]]$estimate), model = "Neonatal death", country = "Malawi")
per.mal    <- cbind(as.data.table(result1[["per.mal"]]$estimate), model = "Pernatal death", country = "Malawi")

still.tan  <- cbind(as.data.table(result1[["still.tan"]]$estimate), model = "All stillbirth", country = "Tanzania")
intra.tan  <- cbind(as.data.table(result1[["intra.tan"]]$estimate), model = "Intra stillbirth", country = "Tanzania")
ante.tan   <- cbind(as.data.table(result1[["ante.tan"]]$estimate), model = "Antermortum", country = "Tanzania")
neo.tan    <- cbind(as.data.table(result1[["neo.tan"]]$estimate), model = "Neonatal death", country = "Tanzania")
per.tan    <- cbind(as.data.table(result1[["per.tan"]]$estimate), model = "Pernatal death", country = "Tanzania")

still.ug  <- cbind(as.data.table(result1[["still.ug"]]$estimate), model = "All stillbirth", country = "Uganda")
intra.ug  <- cbind(as.data.table(result1[["intra.ug"]]$estimate), model = "Intra stillbirth", country = "Uganda")
ante.ug   <- cbind(as.data.table(result1[["ante.ug"]]$estimate), model = "Antermortum", country = "Uganda")
neo.ug    <- cbind(as.data.table(result1[["neo.ug"]]$estimate), model = "Neonatal death", country = "Uganda")
per.ug    <- cbind(as.data.table(result1[["per.ug"]]$estimate), model = "Pernatal death", country = "Uganda")

all <- rbind(still.ben,intra.ben,ante.ben, neo.ben, per.ben,
             still.mal,intra.mal,ante.mal, neo.mal, per.mal,
             still.tan,intra.tan,ante.tan, neo.tan, per.tan,
             still.ug,intra.ug,ante.ug, neo.ug, per.ug)

###### CREATE META-ANALYSES RESULTS
sum <- all[ ( perc == "99th")]

## Create table 4: descriptive PCA analyses
model_st <- c("All stillbirth", "Intra stillbirth", "Antermortum", "Neonatal death", "Pernatal death")
cen <- c("min","25th", "50th", "75th")

# Pooled Estimate
model <- data.frame()  
for (i in model_st){ 
  for (j in cen){  
    
    #i == "All stillbirth"
    #j == "min"  
    
    # Pooled Estimate
    cohort <- sum[model == paste0(i) & cen == paste0(j),.(cen, perc, OR, CIlow, CIhigh, model, country)]
    
    meta_mmt <- metagen(log(OR), lower = log(CIlow), upper = log(CIhigh),
                        studlab = country, data=cohort, sm = "OR")
    
    # extract
    res <-data.table( beta=meta_mmt$TE.random, se=meta_mmt$seTE.random, i2=meta_mmt$I2)
    res[, model := i]
    res[, cen := j]
    model <- as.data.table(rbind(model, res))
  }
}

### make table
model$low<-model$beta-qnorm(1-0.05/2)*model$se
model$high<-model$beta+qnorm(1-0.05/2)*model$se
model$or<-exp(model$beta)
model$or_low<-exp(model$low)
model$or_high<-exp(model$high)


model[, or := formatC(or,digits=2, format="f")]
model[, or_low := formatC(or_low,digits=2, format="f")]
model[, or_high := formatC(or_high, digits=2, format="f")]
model[,c('hr_all'):=paste0(or," (", or_low, "; ",or_high,")")]

# add I2
model[, i2 := formatC(as.numeric(i2),digits=2, format="f")]
save(model, file = "04_results/tmean_lag01_9x9.Rda")

#########################################################################################################
##### Lag 0-2
######################################################################################################### 


### 1. Effect estimation (refer to the function "casecrs" in "CHAMNHA_functions.R") ####
result1 <- list()
for (i in 1:length(dlist)){
  # i = 1
  dat <- dlist[[i]]
  dat <- as.data.frame(dat)
  result1[[i]] <- casecrs_0(status="status", id="alert_id",  lag=2, varper=c(10,75,90), lagnk=2, cen=list(min=TRUE, pct=c(25,50,75)), estpct=c(1,2.5,5,10,90,95,97.5,99), data=dat)
  print(i)
}

names(result1) <- name_outcome


#### Put together the main results
still.ben  <- cbind(as.data.table(result1[["still.ben"]]$estimate), model = "All stillbirth", country = "Benin")
intra.ben  <- cbind(as.data.table(result1[["intra.ben"]]$estimate), model = "Intra stillbirth", country = "Benin")
ante.ben   <- cbind(as.data.table(result1[["ante.ben"]]$estimate), model = "Antermortum", country = "Benin")
neo.ben    <- cbind(as.data.table(result1[["neo.ben"]]$estimate), model = "Neonatal death", country = "Benin")
per.ben    <- cbind(as.data.table(result1[["per.ben"]]$estimate), model = "Pernatal death", country = "Benin")

still.mal  <- cbind(as.data.table(result1[["still.mal"]]$estimate), model = "All stillbirth", country = "Malawi")
intra.mal  <- cbind(as.data.table(result1[["intra.mal"]]$estimate), model = "Intra stillbirth", country = "Malawi")
ante.mal   <- cbind(as.data.table(result1[["ante.mal"]]$estimate), model = "Antermortum", country = "Malawi")
neo.mal    <- cbind(as.data.table(result1[["neo.mal"]]$estimate), model = "Neonatal death", country = "Malawi")
per.mal    <- cbind(as.data.table(result1[["per.mal"]]$estimate), model = "Pernatal death", country = "Malawi")

still.tan  <- cbind(as.data.table(result1[["still.tan"]]$estimate), model = "All stillbirth", country = "Tanzania")
intra.tan  <- cbind(as.data.table(result1[["intra.tan"]]$estimate), model = "Intra stillbirth", country = "Tanzania")
ante.tan   <- cbind(as.data.table(result1[["ante.tan"]]$estimate), model = "Antermortum", country = "Tanzania")
neo.tan    <- cbind(as.data.table(result1[["neo.tan"]]$estimate), model = "Neonatal death", country = "Tanzania")
per.tan    <- cbind(as.data.table(result1[["per.tan"]]$estimate), model = "Pernatal death", country = "Tanzania")

still.ug  <- cbind(as.data.table(result1[["still.ug"]]$estimate), model = "All stillbirth", country = "Uganda")
intra.ug  <- cbind(as.data.table(result1[["intra.ug"]]$estimate), model = "Intra stillbirth", country = "Uganda")
ante.ug   <- cbind(as.data.table(result1[["ante.ug"]]$estimate), model = "Antermortum", country = "Uganda")
neo.ug    <- cbind(as.data.table(result1[["neo.ug"]]$estimate), model = "Neonatal death", country = "Uganda")
per.ug    <- cbind(as.data.table(result1[["per.ug"]]$estimate), model = "Pernatal death", country = "Uganda")

all <- rbind(still.ben,intra.ben,ante.ben, neo.ben, per.ben,
             still.mal,intra.mal,ante.mal, neo.mal, per.mal,
             still.tan,intra.tan,ante.tan, neo.tan, per.tan,
             still.ug,intra.ug,ante.ug, neo.ug, per.ug)

###### CREATE META-ANALYSES RESULTS
sum <- all[ ( perc == "99th")]

## Create table 4: descriptive PCA analyses
model_st <- c("All stillbirth", "Intra stillbirth", "Antermortum", "Neonatal death", "Pernatal death")
cen <- c("min","25th", "50th", "75th")

# Pooled Estimate
model <- data.frame()  
for (i in model_st){ 
  for (j in cen){  
    
    #i == "All stillbirth"
    #j == "min"  
    
    # Pooled Estimate
    cohort <- sum[model == paste0(i) & cen == paste0(j),.(cen, perc, OR, CIlow, CIhigh, model, country)]
    
    meta_mmt <- metagen(log(OR), lower = log(CIlow), upper = log(CIhigh),
                        studlab = country, data=cohort, sm = "OR")
    
    # extract
    res <-data.table( beta=meta_mmt$TE.random, se=meta_mmt$seTE.random, i2=meta_mmt$I2)
    res[, model := i]
    res[, cen := j]
    model <- as.data.table(rbind(model, res))
  }
}
### make table
model$low<-model$beta-qnorm(1-0.05/2)*model$se
model$high<-model$beta+qnorm(1-0.05/2)*model$se
model$or<-exp(model$beta)
model$or_low<-exp(model$low)
model$or_high<-exp(model$high)


model[, or := formatC(or,digits=2, format="f")]
model[, or_low := formatC(or_low,digits=2, format="f")]
model[, or_high := formatC(or_high, digits=2, format="f")]
model[,c('hr_all'):=paste0(or," (", or_low, "; ",or_high,")")]

# add I2
model[, i2 := formatC(as.numeric(i2),digits=2, format="f")]
save(model, file = "04_results/tmean_lag02_9x9.Rda")





