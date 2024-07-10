##############################################################################
#
# Contact: Jeroen de Bont (jeroen.de.bont@ki.se)
#
##############################################################################
#
# Script:  SUMMARIZING FILE WITH ALL THE ANALYSIS
#
##############################################################################

#rm(list=ls())
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
library(ggplot2)
library(openxlsx)
library(meta)


#########################################################################################################
#### Creating table S2
#########################################################################################################

####### Country specific  ####
# load all previous results from the other estiamtes and combine them
load("results_country_still_9x9.Rda")
tmean <- all
load("results_country_max_9x9.Rda")
tmax <- all
load("results_country_minimum_9x9.Rda")
tmin <- all

#### create identifier
tmean[, exp := "tmean"]
tmax[, exp := "tmax"]
tmin[, exp := "tmin"]

#### merge
all <- rbind(tmean, tmax, tmin)
all
# meta analyses
library(meta)

###### CREATE META-ANALYSES RESULTS - 75 - 99th ####
#sum <- all[(cen == "min" ) & ( perc == "95th")]
sum <- all[(cen == "75th" ) & ( perc == "99th")] #prueba
sum[, cen_min_perc := formatC(cen_min_perc,digits=1, format="f")]
sum[, cen_min := paste0(cen_min, " (",cen_min_perc,")" )]

## Create table 4: descriptive PCA analyses
model_st <- c("All stillbirth", "Intra stillbirth", "Antermortum", "Neonatal death", "Pernatal death")
exposure <- c("tmean", "tmax", "tmin")

# Pooled Estimate
model <- data.frame()  
for (i in model_st){  
  #i = "All stillbirth"
  # Pooled Estimate
  cohort <- sum[model == paste0(i),.(cen, perc, OR, CIlow, CIhigh, model, country, exp)]
  
  for (j in exposure){  
    #j = "tmean"  
    cohort_temp <- cohort[exp == paste0(j)]
    
    meta_mmt <- metagen(log(OR), lower = log(CIlow), upper = log(CIhigh),
                        studlab = country, data=cohort_temp, sm = "OR")
    
    # extract
    res <-data.table( beta=meta_mmt$TE.random, se=meta_mmt$seTE.random, i2=meta_mmt$I2)
    res[, model := i]
    res[, exp := j]
    model <- as.data.table(rbind(model, res))
  }
}

### make table
model$low<-model$beta-qnorm(1-0.05/2)*model$se
model$high<-model$beta+qnorm(1-0.05/2)*model$se
model$or<-exp(model$beta)
model$or_low<-exp(model$low)
model$or_high<-exp(model$high)

### merge to create figure
model[, country := "Meta-analysis"]
model[, cen_75 := ""]
db1 <- model[,.(country, model, exp, or, or_low, or_high,cen_75)]
db2 <- sum[,.(country, model, exp, OR, CIlow, CIhigh, cen_75)]
colnames(db2) <- c("country", "model", "exp", "or", "or_low", "or_high", "cen_75")
db <- rbind(db1, db2)

db[, or := formatC(or,digits=2, format="f")]
db[, or_low := formatC(or_low,digits=2, format="f")]
db[, or_high := formatC(or_high, digits=2, format="f")]
db[, cen_75 := formatC(as.numeric(cen_75), digits=1, format="f")]
db[,c('hr_all'):=paste0(or," (", or_low, "; ",or_high,")")]

# cahnge labels
db$exp<- ordered(db$exp, levels = c( "tmean", "tmax","tmin"))
db$country<- ordered(db$country, levels = c("Benin","Malawi" ,"Tanzania", "Uganda",  "Meta-analysis"))

table <- dcast(db, country + exp ~ model, value.var="hr_all") # MAKE NEW TABLE
table2 <- dcast(db, country + exp ~ model, value.var="cen_75") # MAKE NEW TABLE
setDT(table)
colnames(table2) <- c("country", "exp", "mmt_still", "mmt_ante", "mmt_intra", "mmt_neo", "mmt_per")

tab <- merge(table, table2, by = c("country", "exp"))
tab <- tab[,c(1,2,8,3,9,4,10,5,11,6,12,7)]
tab <- tab[order(exp)]
tab

library(openxlsx)
write.xlsx(tab,file="04_results/results_sens_lag06_75-99.xlsx") 

#########################################################################################################
#### Creating table S3
#########################################################################################################

####### meta analyses the results
load("results_country_still_9x9.Rda")
tmean <- all
load("results_country_max_9x9.Rda")
tmax <- all
load("results_country_minimum_9x9.Rda")
tmin <- all

#### create identifier
tmean[, exp := "tmean"]
tmax[, exp := "tmax"]
tmin[, exp := "tmin"]

#### merge
all <- rbind(tmean, tmax, tmin)
all

# meta analyses
library(meta)

###### CREATE META-ANALYSES RESULTS
sum <- all[(cen == "50th" ) & ( perc == "99th")]

## Create table 4: descriptive PCA analyses
model_st <- c("All stillbirth", "Intra stillbirth", "Antermortum", "Neonatal death", "Pernatal death")
exposure <- c("tmean", "tmax", "tmin")
# Pooled Estimate
model <- data.frame()  
for (i in model_st){  
  #i = "All stillbirth"
  # Pooled Estimate
  cohort <- sum[model == paste0(i),.(cen, perc, OR, CIlow, CIhigh, model, country, exp)]
  
  for (j in exposure){  
    #j = "tmean"  
    cohort_temp <- cohort[exp == paste0(j)]
    
    meta_mmt <- metagen(log(OR), lower = log(CIlow), upper = log(CIhigh),
                        studlab = country, data=cohort_temp, sm = "OR")
    
    # extract
    res <-data.table( beta=meta_mmt$TE.random, se=meta_mmt$seTE.random, i2=meta_mmt$I2)
    res[, model := i]
    res[, exp := j]
    model <- as.data.table(rbind(model, res))
  }
}

### make table
model$low<-model$beta-qnorm(1-0.05/2)*model$se
model$high<-model$beta+qnorm(1-0.05/2)*model$se
model$or<-exp(model$beta)
model$or_low<-exp(model$low)
model$or_high<-exp(model$high)

### merge to create figure
model[, country := "Meta-analysis"]
model[, cen_50 := ""]
db1 <- model[,.(country, model, exp, or, or_low, or_high, cen_50)]
db2 <- sum[,.(country, model, exp, OR, CIlow, CIhigh, cen_50)]
colnames(db2) <- c("country", "model", "exp", "or", "or_low", "or_high", "cen_50")
db <- rbind(db1, db2)

db[, or := formatC(or,digits=2, format="f")]
db[, or_low := formatC(or_low,digits=2, format="f")]
db[, or_high := formatC(or_high, digits=2, format="f")]
db[, cen_50 := formatC(as.numeric(cen_50), digits=1, format="f")]
db[,c('hr_all'):=paste0(or," (", or_low, "; ",or_high,")")]

# cahnge labels
db$exp<- ordered(db$exp, levels = c( "tmean", "tmax","tmin"))
db$country<- ordered(db$country, levels = c("Benin","Malawi" ,"Tanzania", "Uganda",  "Meta-analysis"))

table <- dcast(db, country + exp ~ model, value.var="hr_all") # MAKE NEW TABLE
table2 <- dcast(db, country + exp ~ model, value.var="cen_50") # MAKE NEW TABLE
setDT(table)
colnames(table2) <- c("country", "exp", "50th_still", "50th_ante", "50th_intra", "50th_neo", "mmt_per")

tab <- merge(table, table2, by = c("country", "exp"))
tab <- tab[,c(1,2,8,3,9,4,10,5,11,6,12,7)]
tab <- tab[order(exp)]
tab
library(openxlsx)
write.xlsx(tab,file="04_results/results_sens_lag06_p50.xlsx") 

#########################################################################################################
#### Heat effects during 6 consecutive hottest months (Figure S3)
#########################################################################################################

#### summarize all the effect estimates
load("04_results/tmean_lag06_warm_season_9x9.Rda")
tmean <- model
load("04_results/tmax_lag06_warm_season_9x9.Rda")
tmax <- model
load("04_results/min_lag06_warm_season_9x9.Rda")
tmin <- model

#### create identifier
tmean[, exp := "tmean"]
tmax[, exp := "tmax"]
tmin[, exp := "tmin"]

#### merge
all <- rbind(tmean, tmax, tmin)
all$exp<- ordered(all$exp, levels = c("tmean", "tmax" , "tmin"))

table <- dcast(all, exp + cen ~ model, value.var="hr_all") # MAKE NEW TABLE
table
# all
write.xlsx(table,file="04_results/results_lag06_all_summer.xlsx") 


#### heat effects
min <- all
min[, or := as.numeric(or)]
min[, or_low := as.numeric(or_low)]
min[, or_high := as.numeric(or_high)]
min$model<- ordered(min$model, levels = c("Pernatal death","Neonatal death", "Antermortum", "Intra stillbirth", "All stillbirth"))
min$model<- factor(min$model, labels = c("Perinatal death","Neonatal death" , "Macerated stillbirth","Fresh stillbirth",  "Stillbirth"))
min$meas<- ordered(min$exp, labels = c("Mean temperature","Max temperature", "Min temperature"))

## create figure with odds ratios
db_sub <- min

g1 <- ggplot(db_sub, aes(x=cen, y=or, ymin=or_low, ymax=or_high, col=cen, fill=cen)) + 
  geom_linerange(size=1,position=position_dodge(width = 0.5)) +
  geom_hline(yintercept=1, lty=2) +
  geom_point(size=1.5, shape=21, colour="white", stroke = 0.5,position=position_dodge(width = 0.5)) +
  scale_fill_grey( name = "Centered")+
  scale_color_grey( name = "Centered") +
  coord_flip(ylim=c(0, 5)) + theme_minimal() +
  facet_grid(model~meas) +
  ylab("OR (95%CI)") + 
  guides(fill = guide_legend(reverse = TRUE), col = guide_legend(reverse = TRUE)) + theme_bw() +
  theme(axis.title.y=element_blank())

g1
ggsave("05_figures/mmt_centered_summer.tiff", dpi = 300)

##### compare summer vs whole year
load("04_results/tmean_lag06_9x9.Rda")
tmean <- model

load("04_results/tmean_lag06_warm_season_9x9.Rda")
tmean_sum <- model

tmean[, time := "Whole year"]
tmean_sum[, time := "Hottest months"]

all <- rbind(tmean, tmean_sum)
all <- all[ cen == "75th" ]

## all
all[, or := as.numeric(or)]
all[, or_low := as.numeric(or_low)]
all[, or_high := as.numeric(or_high)]

#### heat effects
all$model<- ordered(all$model, levels = c("Pernatal death","Neonatal death", "Intra stillbirth", "Antermortum",  "All stillbirth"))
all$model<- factor(all$model, labels = c("Perinatal death","Neonatal death" , "Intrapartum stillbirth", "Antepartum stillbirth",  "Stillbirth"))

db_sub <- all

g1 <- ggplot(db_sub, aes(x=model, y=or, ymin=or_low, ymax=or_high, col=time, fill=time)) + 
  geom_linerange(size=1,position=position_dodge(width = 0.5)) +
  geom_hline(yintercept=1, lty=2) +
  geom_point(size=1.5, shape=21, colour="white", stroke = 0.5,position=position_dodge(width = 0.5)) +
  scale_fill_grey(name = "")+
  scale_color_grey(name = "") +
  coord_flip(ylim=c(0, 5)) + theme_minimal() +
  ylab("OR (95%CI)") + 
  guides(fill = guide_legend(reverse = TRUE), col = guide_legend(reverse = TRUE)) + theme_bw() +
  theme(axis.title.y=element_blank()) +
  theme(text=element_text(size=14,  family="sans"))

g1
ggsave("05_figures/summer months.tiff",dpi = 300)

#########################################################################################################
#### Lag comparison (Figure S4)
#########################################################################################################

## temperature main tables
load("04_results/tmean_lag0_9x9.Rda")
lag0 <- model
load("04_results/tmean_lag01_9x9.Rda")
lag01 <- model
load("04_results/tmean_lag02_9x9.Rda")
lag02 <- model
load("04_results/tmean_lag06_9x9.Rda")
lag06 <- model

# add values
lag0[, lag := "lag0"]
lag01[, lag := "lag01"]
lag02[, lag := "lag02"]
lag06[, lag := "lag06"]

#merge
temp <- rbind(lag0, lag01, lag02, lag06)
temp[, meas := "Mean temperature"]

## Max temperature main tables
load("04_results/tmax_lag0_9x9.Rda")
lag0 <- model
load("04_results/tmax_lag01_9x9.Rda")
lag01 <- model
load("04_results/tmax_lag02_9x9.Rda")
lag02 <- model
load("04_results/tmax_lag06_9x9.Rda")
lag06 <- model

# add values
lag0[, lag := "lag0"]
lag01[, lag := "lag01"]
lag02[, lag := "lag02"]
lag06[, lag := "lag06"]

temp_max <- rbind(lag0, lag01, lag02, lag06)
temp_max[, meas := "Max temperature"]

## Min temperature main tables
load("04_results/tmin_lag0_9x9.Rda")
lag0 <- model
load("04_results/tmin_lag01_9x9.Rda")
lag01 <- model
load("04_results/tmin_lag02_9x9.Rda")
lag02 <- model
load("04_results/tmin_lag06_9x9.Rda")
lag06 <- model

# add values
lag0[, lag := "lag0"]
lag01[, lag := "lag01"]
lag02[, lag := "lag02"]
lag06[, lag := "lag06"]

temp_min <- rbind(lag0, lag01, lag02, lag06)
temp_min[, meas := "Min temperature"]

all <-  rbind(temp, temp_max, temp_min)
all

#### prepare for figure
min <- all
min[, or := as.numeric(or)]
min[, or_low := as.numeric(or_low)]
min[, or_high := as.numeric(or_high)]
min$lag<- ordered(min$lag, levels = c("lag06" ,"lag02", "lag01",  "lag0"))
min$model<- ordered(min$model, levels = c("Pernatal death","Neonatal death", "Intra stillbirth", "Antermortum",  "All stillbirth"))
min$model<- factor(min$model, labels = c("Perinatal death","Neonatal death" , "Intrapartum stillbirth", "Antepartum stillbirth",  "Stillbirth"))
min$meas<- ordered(min$meas, levels = c("Mean temperature","Max temperature", "Min temperature"))

## create figure with odds ratios
db_sub <- min[cen == "75th" & meas == "Mean temperature" ]

g1 <- ggplot(db_sub, aes(x=model, y=or, ymin=or_low, ymax=or_high, col=lag, fill=lag)) + 
  geom_linerange(size=1,position=position_dodge(width = 0.5)) +
  geom_hline(yintercept=1, lty=2) +
  geom_point(size=1.5, shape=21, colour="white", stroke = 0.5,position=position_dodge(width = 0.5)) +
  scale_fill_grey( name = "Lags")+
  scale_color_grey( name = "Lags") +
  #scale_x_discrete(name="") + 
  #scale_y_continuous(lim = c(0.5,5), breaks = seq(0.5,5,by = 0.5)) +
  coord_flip(ylim=c(0, 4)) + theme_minimal() +
  #facet_grid(~meas) +
   ylab("OR (95%CI)") + 
  guides(fill = guide_legend(reverse = TRUE), col = guide_legend(reverse = TRUE)) + theme_bw() +
  theme(axis.title.y=element_blank()) +
  theme(text=element_text(size=14,  family="sans"))

g1
ggsave("05_figures/mmt_lag.tiff",dpi = 300)

#########################################################################################################
#### Sensitivity analyses by number knots and their location (Figure S5)
#########################################################################################################

### main effect
#### summarize all the effect estuimates
load("04_results/tmean_lag06_9x9.Rda")
main <- model[cen == "75th"]
load("04_results/sens_1k_50.Rda")
kn1_50 <- model
load("04_results/sens_1k_75.Rda")
kn1_75 <- model
load("04_results/sens_5k_25_75.Rda")
kn2_25_75 <- model
load("04_results/sens_5k_50_75.Rda")
kn2_50_75 <- model
load("04_results/sens_5k_50_90.Rda")
kn2_50_90 <- model
load("04_results/sens_5k_75_90.Rda")
kn2_75_90 <- model
load("04_results/sens_5k_25_50_75.Rda")
kn3_25_50_75 <- model
load("04_results/sens_5k_10_50_90.Rda")
kn3_10_50_90 <- model

# make an identifier
main[, sens := "Main model"]
main[, cen := NULL]
kn1_50[, sens := "1 knot 50th"]
kn1_75[, sens := "1 knot 75th"]
kn2_25_75[, sens := "2 knot 25 and 50th"]
kn2_50_75[, sens := "2 knot 50 and 75th"]
kn2_50_90[, sens := "2 knot 50 and 90th"]
kn2_75_90[, sens := "2 knot 75 and 90th"]
kn3_25_50_75[, sens := "3 knot 25, 50 and 75th"]
kn3_10_50_90[, sens := "3 knot 10, 50 and 90th"]

#### merge
min <- rbind(main, kn1_50, kn1_75, kn2_25_75, kn2_50_75, kn2_50_90, kn2_75_90, kn3_25_50_75, kn3_10_50_90)

min[, or := as.numeric(or)]
min[, or_low := as.numeric(or_low)]
min[, or_high := as.numeric(or_high)]

#### heat effects
min$model<- ordered(min$model, levels = c("Pernatal death","Neonatal death", "Intra stillbirth", "Antermortum",  "All stillbirth"))
min$model<- factor(min$model, labels = c("Perinatal death","Neonatal death", "Intrapartum stillbirth", "Antepartum stillbirth",  "Stillbirth"))

db_sub <- min

g1 <- ggplot(db_sub, aes(x=model, y=or, ymin=or_low, ymax=or_high, col=sens, fill=sens)) + 
  geom_linerange(size=1,position=position_dodge(width = 0.5)) +
  geom_hline(yintercept=1, lty=2) +
  geom_point(size=1.5, shape=21, colour="white", stroke = 0.5,position=position_dodge(width = 0.5)) +
  scale_fill_grey(name = "")+
  scale_color_grey(name = "") +
  #scale_x_discrete(name="") + 
  #scale_y_continuous(lim = c(0.5,5), breaks = seq(0.5,5,by = 0.5)) +
  coord_flip(ylim=c(0, 4)) + theme_minimal() +
  # facet_grid(~exp) +
  ylab("OR (95%CI)") + 
  guides(fill = guide_legend(reverse = TRUE), col = guide_legend(reverse = TRUE)) + theme_bw() +
  theme(axis.title.y=element_blank())  +
  theme(text=element_text(size=14,  family="sans"))

g1
ggsave("05_figures/sens_splines.tiff",dpi = 300)

#########################################################################################################
#### Results with and without referral births (Figure S6)
#########################################################################################################

#### Results without refferral ####
### main effect
load("04_results/tmean_lag06_9x9.Rda")
main <- model

### without referral
load("04_results/tmean_lag06_9x9_withouth referral.Rda")
with <- model

#### create identifier
main[, exp := "Main model"]
with[, exp := "Without referral births"]


#### merge
all <- rbind(main,with)
all <- all[ cen == "75th"] 

## all
min <-  all
min[, or := as.numeric(or)]
min[, or_low := as.numeric(or_low)]
min[, or_high := as.numeric(or_high)]

#### heat effects
min$model<- ordered(min$model, levels = c("Pernatal death", "Intra stillbirth", "Antermortum",  "All stillbirth"))
min$model<- factor(min$model, labels = c("Perinatal death", "Intrapartum stillbirth", "Antepartum stillbirth",  "Stillbirth"))

db_sub <- min

g1 <- ggplot(db_sub, aes(x=model, y=or, ymin=or_low, ymax=or_high, col=exp, fill=exp)) + 
  geom_linerange(size=1,position=position_dodge(width = 0.5)) +
  geom_hline(yintercept=1, lty=2) +
  geom_point(size=1.5, shape=21, colour="white", stroke = 0.5,position=position_dodge(width = 0.5)) +
  scale_fill_grey(name = "")+
  scale_color_grey(name = "") +
  #scale_x_discrete(name="") + 
  #scale_y_continuous(lim = c(0.5,5), breaks = seq(0.5,5,by = 0.5)) +
  coord_flip(ylim=c(0, 5)) + theme_minimal() +
  # facet_grid(~exp) +
  ylab("OR (95%CI)") + 
  guides(fill = guide_legend(reverse = TRUE), col = guide_legend(reverse = TRUE)) + theme_bw() +
  theme(axis.title.y=element_blank()) +
  theme(text=element_text(size=14,  family="sans"))


g1

ggsave("05_figures/mean_without referral.tiff",dpi = 300)


#########################################################################################################
#### Adjusting for relative humidity (Figure S7)
#########################################################################################################

#### Relative humidity ####
### main effect
#### summarize all the effect estuimates
load("04_results/tmean_lag06_9x9.Rda")
tmean <- model
load("04_results/tmax_lag06_9x9.Rda")
tmax <- model
load("04_results/tmin_lag06_9x9.Rda")
tmin <- model

#### create identifier
tmean[, exp := "tmean"]
tmax[, exp := "tmax"]
tmin[, exp := "tmin"]

#### merge
all <- rbind(tmean, tmax, tmin)
all[, meas := "Not adj. RH" ]
all <- all[ cen == "75th"]

### the dew adjusted models
load("04_results/sens_tmean_lag06_9x9.Rda")
tmean_h <- model
load("04_results/sens_tmax_lag06_9x9.Rda")
tmax_h <- model
load("04_results/sens_tmin_lag06_9x9.Rda")
tmin_h <- model

#### create identifier
tmean_h[, exp := "tmean"]
tmax_h[, exp := "tmax"]
tmin_h[, exp := "tmin"]

#### merge
dew <- rbind(tmean_h, tmax_h, tmin_h)
dew[, meas := "adj. RH" ]

all <- all[,.(beta, se, i2, model, low, high, or, or_low, or_high, hr_all, exp, meas)]

## all
min <-  rbind(all, dew)
min[, or := as.numeric(or)]
min[, or_low := as.numeric(or_low)]
min[, or_high := as.numeric(or_high)]

#### heat effects
min$model<- ordered(min$model, levels = c("Pernatal death","Neonatal death", "Intra stillbirth", "Antermortum",  "All stillbirth"))
min$model<- factor(min$model, labels = c("Perinatal death","Neonatal death", "Intrapartum stillbirth", "Antepartum stillbirth",  "Stillbirth"))

min$exp<- factor(min$exp, labels = c("Max temperature","Mean temperature" , "Min temperature"))
min$exp<- ordered(min$exp, levels = c("Mean temperature","Max temperature", "Min temperature"))

db_sub <- min[model != "Neonatal death"]
db_sub <- db_sub[exp == "Mean temperature"]

g1 <- ggplot(db_sub, aes(x=model, y=or, ymin=or_low, ymax=or_high, col=meas, fill=meas)) + 
  geom_linerange(size=1,position=position_dodge(width = 0.5)) +
  geom_hline(yintercept=1, lty=2) +
  geom_point(size=1.5, shape=21, colour="white", stroke = 0.5,position=position_dodge(width = 0.5)) +
  scale_fill_grey(name = "")+
  scale_color_grey(name = "") +
  #scale_x_discrete(name="") + 
  #scale_y_continuous(lim = c(0.5,5), breaks = seq(0.5,5,by = 0.5)) +
  coord_flip(ylim=c(0,4)) + theme_minimal() +
  # facet_grid(~exp) +
  ylab("OR (95%CI)") + 
  guides(fill = guide_legend(reverse = TRUE), col = guide_legend(reverse = TRUE)) + theme_bw() +
  theme(axis.title.y=element_blank()) +
  theme(text=element_text(size=14,  family="sans"))

g1
ggsave("05_figures/mean_adj_hum.tiff",dpi = 300)


#########################################################################################################
#### MAIN ANALYSES USING THE 9X9KM GRIDDED TEMPERATURE VS 28X28KM GRIDDED TEMPERATURE (FIGURE S8)
#########################################################################################################

# Set directory
setwd("/Volumes/projects$/EnvEpi-Projekt/CHAMNHA/4_paper_alert")
sink("log_file_results.txt")

#### summarize all the effect estimates from mean, max and min temperature results + 28x28km results
# 9x9
load("04_results/tmean_lag06_9x9.Rda")
tmean_9x9 <- model
load("04_results/tmax_lag06_9x9.Rda")
tmax_9x9 <- model
load("04_results/tmin_lag06_9x9.Rda")
tmin_9x9 <- model

# 28x28
load("04_results/tmean_lag06_28x28.Rda")
tmean_28x28 <- model
load("04_results/tmax_lag06_28x28.Rda")
tmax_28x28 <- model
load("04_results/tmin_lag06_28x28.Rda")
tmin_28x28 <- model

#### create identifier
tmean_9x9[, exp := "tmean"]
tmax_9x9[, exp := "tmax"]
tmin_9x9[, exp := "tmin"]

tmean_28x28[, exp := "tmean"]
tmax_28x28[, exp := "tmax"]
tmin_28x28[, exp := "tmin"]

#### merge
all_9 <- rbind(tmean_9x9, tmax_9x9, tmin_9x9)
all_9$exp<- ordered(all_9$exp, levels = c("tmean", "tmax" , "tmin"))

all_28<- rbind(tmean_28x28, tmax_28x28, tmin_28x28)
all_28$exp<- ordered(all_28$exp, levels = c("tmean", "tmax" , "tmin"))

# combine
all_9[,resolution := "9x9km"]
all_28[,resolution := "28x28km"]

all <- rbind(all_9, all_28)

# create table
table <- dcast(all, exp + cen + resolution~ model, value.var="hr_all") # MAKE NEW TABLE
table

# all
write.xlsx(table,file="04_results/results_lag06_9x9vs28x28.xlsx") 

#only 75th percentile as centered value
tab <- table[cen == "75th"]

write.xlsx(tab,file="04_results/results_lag06_75th_9x9vs28x28.xlsx") 

### Create figure

min <- all[cen == "75th"]
min[, or := as.numeric(or)]
min[, or_low := as.numeric(or_low)]
min[, or_high := as.numeric(or_high)]
min$model<- ordered(min$model, levels = c("Pernatal death","Neonatal death", "Antermortum", "Intra stillbirth", "All stillbirth"))
min$model<- factor(min$model, labels = c("Perinatal death","Neonatal death" , "Macerated stillbirth","Fresh stillbirth",  "Stillbirth"))
min$meas<- ordered(min$exp, labels = c("Mean temperature","Max temperature", "Min temperature"))

## create figure with odds ratios
# remove neonatal deaths
db_sub <- min[model != "Neonatal death"]
db_sub <- db_sub[meas == "Mean temperature"]

g1 <- ggplot(db_sub, aes(x=model, y=or, ymin=or_low, ymax=or_high, col=resolution, fill=resolution)) + 
  geom_linerange(size=1,position=position_dodge(width = 0.5)) +
  geom_hline(yintercept=1, lty=2) +
  geom_point(size=1.5, shape=21, colour="white", stroke = 0.5,position=position_dodge(width = 0.5)) +
  scale_fill_grey( name = "Resolution:")+
  scale_color_grey( name = "Resolution:") +
  coord_flip(ylim=c(0, 5)) + theme_minimal() +
  #facet_grid(~model) +
  ylab("OR (95%CI)") + 
  guides(fill = guide_legend(reverse = TRUE), col = guide_legend(reverse = TRUE)) + theme_bw() +
  theme(axis.title.y=element_blank()) +
  theme(text=element_text(size=14,  family="sans"))


g1
ggsave("05_figures/9x9km vs 28x28km.tiff",dpi = 300)

#########################################################################################################
#### Long-term effects by trimester-specific exposure (Figure S9)
#########################################################################################################
 
# 9x9
load("04_results/long_mean_9x9.Rda")
tmean_9x9 <- model
load("04_results/long_min_9x9.Rda")
tmax_9x9 <- model
load("04_results/long_max_9x9.Rda")
tmin_9x9 <- model


#### merge
all<- rbind(tmean_9x9, tmax_9x9, tmin_9x9)
all$out<- ordered(all$out, levels = c("stil","ante", "intr", "peri"))
all$temp<- ordered(all$temp, levels = c("Mean","Min", "Max"))

table <- dcast(all,  model + out + country ~temp , value.var="hr_all") # MAKE NEW TABLE
table

# all
write.xlsx(table,file="04_results/result_long_all.xlsx") 

#### heat effects: focus on 9x9lm
min <- all[temp == "Mean"]
min[, or := as.numeric(or)]
min[, or_low := as.numeric(or_low)]
min[, or_high := as.numeric(or_high)]
min$out<- ordered(min$out, levels = c( "peri", "intr", "ante",   "stil"))
min$out<- factor(min$out, labels = c("Perinatal death", "Intrapartum stillbirths", "Antepartum stillbirths", "Stillbirth"))
min$country<- ordered(min$country, levels = c("Third", "Second", "First"))
#min$temp<- ordered(min$temp, levels = c("Mean","Min", "Max"))


## create figure with odds ratios
# remove neonatal deaths
cox <- min[model == "Cox"]
g1 <- ggplot(cox, aes(x=out, y=or, ymin=or_low, ymax=or_high, col=country, fill=country)) + 
  geom_linerange(size=1,position=position_dodge(width = 0.5)) +
  geom_hline(yintercept=1, lty=2) +
  geom_point(size=1.5, shape=21, colour="white", stroke = 0.5,position=position_dodge(width = 0.5)) +
  scale_fill_grey( name = "Trimester:")+
  scale_color_grey( name = "Trimester:") +
  coord_flip() + theme_minimal() +
  #facet_grid(~temp) +
  ylab("HR (95%CI)") + 
  guides(fill = guide_legend(reverse = TRUE), col = guide_legend(reverse = TRUE)) + theme_bw() +
  theme(axis.title.y=element_blank()) +
  theme(text=element_text(size=14,  family="sans"))

g1
ggsave("05_figures/long_term_cox.tiff",dpi = 300)


#########################################################################################################
#### Figures S10 and S11 can be created using the 03_temperature_stillbirth_meta-regression.Rcode
#### this codes needs to be substituted with maximum and minimum temperature
#########################################################################################################

#########################################################################################################
#### STRATIFIED ANALYSES (Figure 12)
#########################################################################################################

#### summarize all the effect estimates
load("04_results/tmean_lag06_girls_9.Rda")
girl <- model
load("04_results/tmean_lag06_boy_9.Rda")
boy <- model
load("04_results/tmean_lag06_preg_num2+_9.Rda")
pre_num2 <- model
load("04_results/tmean_lag06_preg_num1_9.Rda")
pre_num1 <- model
load("04_results/tmean_lag06_no_preterm_9.Rda")
nopre <- model
load("04_results/tmean_lag06_preterm_9.Rda")
pre <- model
load("04_results/tmean_lag06_mat>35_9.Rda")
age35 <- model
load("04_results/tmean_lag06_mat<35_9.Rda")
agelow_35 <- model
load("04_results/tmean_lag06_nolowbw_9.Rda")
lbw <- model
load("04_results/tmean_lag06_lowbirth_yes_9.Rda")
no_lbw <- model
load("04_results/tmean_lag06_preg_hot_9.Rda")
hot <- model
load("04_results/tmean_lag06_hiv_no_9.Rda")
hiv_no <- model
load("04_results/tmean_lag06_hiv_yes_9.Rda")
hiv_ys <- model
load("04_results/tmean_lag06_obs_no_9.Rda")
obs_no <- model
load("04_results/tmean_lag06_obs_yes_9.Rda")
obs_ys <- model
load("04_results/tmean_lag06_aph_no_9.Rda")
aph_no <- model
load("04_results/tmean_lag06_aph_yes_9.Rda")
aph_ys <- model
load("04_results/tmean_lag06_ht_no_9.Rda")
ht_no <- model
load("04_results/tmean_lag06_ht_yes_9.Rda")
ht_ys <- model

#### create identifier
girl[, eff := "Girl"]
boy[, eff := "Boy"]
pre_num1[, eff := "1st"]
pre_num2[, eff := "2nd or more"]
pre[, eff := "Preterm"]
nopre[, eff := "To term"]
agelow_35[, eff := "<35 years"]
age35[, eff := ">=35 years"]
lbw[, eff := "<2500 gr"]
no_lbw[, eff := ">= 2500 gr"]
hot[, eff := "Hot"]
hiv_no[, eff := "Negative"]
hiv_ys[, eff := "Positive"]
obs_no[, eff := "Normal"]
obs_ys[, eff := "Prolongued "]
aph_no[, eff := "No"]
aph_ys[, eff := "Yes"]
ht_no[, eff := "No"]
ht_ys[, eff := "Yes"]

girl[, eff_mod := "Sex"]
boy[, eff_mod := "Sex"]
pre_num1[, eff_mod := "N° pregnancies"]
pre_num2[, eff_mod := "N° pregnancies"]
pre[, eff_mod := "Preterm"]
nopre[, eff_mod := "Preterm"]
agelow_35[, eff_mod := "Maternal age"]
age35[, eff_mod := "Maternal age"]
lbw[, eff_mod := "Birthweight"]
no_lbw[, eff_mod := "Birthweight"]
hot[, eff_mod := "Season"]
hiv_no[, eff_mod := "HIV status"]
hiv_ys[, eff_mod := "HIV status"]
obs_no[, eff_mod := "Labour problems"]
obs_ys[, eff_mod := "Labour problems"]
aph_no[, eff_mod := "APH"]
aph_ys[, eff_mod := "APH"]
ht_no[, eff_mod := "Hypertensive"]
ht_ys[, eff_mod := "Hypertensive"]

#### merge
all <- rbind(girl, boy, pre_num1, pre_num2,pre, nopre, agelow_35, age35, lbw, no_lbw, hiv_no, hiv_ys, obs_no, obs_ys, aph_no, aph_ys, ht_no, ht_ys)
all$model<- ordered(all$model, levels = c("All stillbirth", "Antermortum",  "Intra stillbirth", "Pernatal death","Neonatal death"))
all$model<- factor(all$model, labels = c("Stillbirth","Antepartum stillbirth","Intrapartum stillbirth", "Perinatal death","Neonatal death"))
all$eff_mod<- ordered(all$eff_mod, levels = c("Season",  "Maternal age", "HIV status", "Hypertensive", "N° pregnancies", "Sex",
                                              "Preterm", "Birthweight", "Labour problems", "APH"))
## create figure with odds ratios
# keep only relevatn
db_sub <- all[model != "Neonatal death"]
#db_sub <- db_sub[model != "Fresh stillbirth"]
#db_sub <- db_sub[model != "Macerated stillbirth"]

db_sub[, or := as.numeric(or)]
db_sub[, or_low := as.numeric(or_low)]
db_sub[, or_high := as.numeric(or_high)]
db_sub
#db_sub <- db_sub[model != "Perinatal death"]
g1 <- ggplot(db_sub, aes(x=eff, y=or, ymin=or_low, ymax=or_high, col=eff_mod, fill=eff_mod)) + 
  geom_linerange(size=1,position=position_dodge(width = 0.5)) +
  geom_hline(yintercept=1, lty=2) +
  geom_point(size=2, shape=21, colour="white", stroke = 0.5,position=position_dodge(width = 0.5)) +
  coord_flip(ylim=c(0, 5)) + theme_minimal() +
  facet_grid(eff_mod~model, scales = "free_y", switch = "y", space = 'free') +
  theme(strip.placement = "outside",
        strip.text.y = element_text(face = "bold", angle=180, vjust = 1))+
  ylab("OR (95%CI)") + scale_fill_grey() + scale_color_grey() +
  guides(fill = guide_legend(reverse = TRUE), col = guide_legend(reverse = TRUE)) + theme_bw() +
  theme(axis.title.y=element_blank(), legend.position="none")+
  theme(text=element_text(size=14,  family="sans"))

g1
ggsave("05_figures/stratified.tiff",width = 30, height = 30, units = "cm", dpi = 300)
table <- dcast(db_sub, eff ~ model, value.var="hr_all") # MAKE NEW TABLE
table
# all
write.xlsx(table,file="04_results/results_lag06_stratified.xlsx") 

sink()
