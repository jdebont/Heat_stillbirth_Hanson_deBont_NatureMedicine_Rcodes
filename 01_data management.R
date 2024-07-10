##############################################################################
#
# Contact: Jeroen de Bont (jeroen.de.bont@ki.se)
#
##############################################################################
#
# Script:  Data management
#
##############################################################################

library(plyr)
library(Hmisc)
library(data.table)
require(compareGroups)
require(dplyr)
require(expss)
library(ggplot2)
library(haven)
library(dlnm)

# setnames
#setwd("/Users/jeroen.de.bont/Desktop/3_paper_alert/")
setwd("/Volumes/projects$/EnvEpi-Projekt/CHAMNHA/4_paper_alert")

#### clean temperature dataset ####
#load the original exposure dataset (28x28km)
tem0 <- readRDS("02_original/ALERT_ID_ERA5_025x025_01012019_31122019.RDS")
tem1 <- readRDS("02_original/ALERT_ID_ERA5_2020_2022.Rds")
tem2 <- readRDS("02_original/ALERT_ID_ERA5_01102022_31122022.Rds")
tem3 <- readRDS("02_original/ALERT_ID_ERA5_01012023_25052023.Rds")
tem4 <- readRDS("02_original/ALERT_ID_ERA5_26052023_05072023.Rds")
tem5 <- readRDS("02_original/ALERT_ID_ERA5_01072023_30092023.Rds")
tem6 <- readRDS("02_original/ALERT_ID_ERA5_025x025_01102023_31122023.RDS")

# new temperature valuues (9x9km)
tem_9x9 <- readRDS("02_original/ID_ALERT_ERA5_temp_01x01_01012019_31122023.Rds")

# check names
names(tem0)
names(tem1)
names(tem2)
names(tem3)
names(tem4)
names(tem5)
names(tem6)
names(tem_9x9)

#change of all variables
setnames(tem0, "2m_dewpoint_temperature", "t2d")
setnames(tem0, "Date", "day")

setnames(tem1, "2m_temperature", "t2m")
setnames(tem1, "2m_dewpoint_temperature", "t2d")
setnames(tem1, "Date", "day")

setnames(tem2, "2m_dewpoint_temperature", "t2d")
setnames(tem2, "Date", "day")

setnames(tem3, "2m_dewpoint_temperature", "t2d")
setnames(tem3, "Date", "day")

setnames(tem4, "2m_dewpoint_temperature", "t2d")
setnames(tem4, "Date", "day")

setnames(tem5, "2m_dewpoint_temperature", "t2d")
setnames(tem5, "Date", "day")

setnames(tem6, "2m_dewpoint_temperature", "t2d")
setnames(tem6, "Date", "day")

setnames(tem_9x9, "2m_dewpoint_temperature", "t2d")
setnames(tem_9x9, "Date", "day")

# keep the same variables we are going to use
tem0 <- tem0[,.(day, alert_id, idera5, t2d)]
tem1 <- tem1[,.(day, alert_id, idera5, t2d)]
tem2 <- tem2[,.(day, alert_id, idera5, t2d)]
tem3 <- tem3[,.(day, alert_id, idera5, t2d)]
tem4 <- tem4[,.(day, alert_id, idera5, t2d)]
tem5 <- tem5[,.(day, alert_id, idera5, t2d)]
tem6 <- tem6[,.(day, alert_id, idera5, t2d)]

# remove duplicates
tem5 <- tem5[day > "2023-07-05"]

# combine 28x28
tem28 <- rbind(tem0, tem1, tem2, tem3, tem4, tem5, tem6)

#check duplicates
tem28[, dup := 1:.N, by = c("day", "alert_id")]
table(tem28$dup) # no duplicates
tem28[, dup := NULL]

# open the other temperature dataset that includes max, min sd, etc
tem0_1 <- readRDS("02_original/ALERT_ID_ERA5_025x025_01012019_31122019.RDS")
tem1_1 <- readRDS("02_original/ALERT_ID_ERA5_MIN_MAX_MEAN_SD_TEMP_2020_2022.Rds")
tem2_1 <- readRDS("02_original/ALERT_ID_ERA5_MIN_MAX_MEAN_SD_TEMP_01112022_31122022.Rds")
tem3_1 <- readRDS("02_original/ALERT_ID_ERA5_01012023_25052023.Rds")
tem4_1 <- readRDS("02_original/ALERT_ID_ERA5_26052023_05072023.Rds")
tem5_1 <- readRDS("02_original/ALERT_ID_ERA5_01072023_30092023.Rds")
tem6_1 <- readRDS("02_original/ALERT_ID_ERA5_025x025_01102023_31122023.RDS")

#keep the relevant variables
tem0_1 <- tem0_1[,.(Date, alert_id, min_temp, max_temp, mean_temp,sd_temp)]
tem1_1 <- tem1_1[,.(Date, alert_id, min_temp, max_temp, mean_temp,sd_temp)]
tem2_1 <- tem2_1[,.(Date, alert_id, min_temp, max_temp, mean_temp,sd_temp)]
tem3_1 <- tem3_1[,.(Date, alert_id, min_temp, max_temp, mean_temp,sd_temp)]
tem4_1 <- tem4_1[,.(Date, alert_id, min_temp, max_temp, mean_temp,sd_temp)]
tem5_1 <- tem5_1[,.(Date, alert_id, min_temp, max_temp, mean_temp,sd_temp)]
tem6_1 <- tem6_1[,.(Date, alert_id, min_temp, max_temp, mean_temp,sd_temp)]

# remove duplicates
tem5_1 <- tem5_1[Date > "2023-07-05"]

setnames(tem0_1, "Date", "day")
setnames(tem1_1, "Date", "day")
setnames(tem2_1, "Date", "day")
setnames(tem3_1, "Date", "day")
setnames(tem4_1, "Date", "day")
setnames(tem5_1, "Date", "day")
setnames(tem6_1, "Date", "day")

tem_mean <- rbind(tem0_1, tem1_1, tem2_1, tem3_1, tem4_1, tem5_1, tem6_1)

# combine exposures of hte 28x28km exposures
all <- merge(tem28, tem_mean, by = c("day", "alert_id"))
setnames(all, "mean_temp", "t2m") # temporary we are going to use this variable

### setnames of the 9x9km combined them
colnames(tem_9x9) <- c("idera5_9", "alert_id", "distance", "day", "t2d_9", "t2m_9" ,"min_temp_9", "max_temp_9", "sd_temp_9", "year")

## merge
all <- merge(all, tem_9x9, by = c("day", "alert_id"), all = TRUE)

## check duplicates:
all[, dup := 1:.N, by = c("day", "alert_id")]
table(all$dup) # no duplicates
all[, dup := NULL]

### estimate percentile levels by hopsital
dt <- all[,.(day, alert_id, t2m_9)]
dt <- dt[order(alert_id, day)]

# Create a new column 'percentile' using the data.table syntax
dt[, perc := as.numeric(ecdf(t2m_9)(t2m_9)), by = alert_id]

# If you want to scale the percentiles from 0 to 100
dt[, perc := perc * 100, by = alert_id]
q<- dt[alert_id ==11] # looks good

dt <- dt[,.(day, alert_id, perc)]

## merge
all <- merge(all, dt, by = c("day", "alert_id"), all = TRUE)

##### Create time-series 2014-2019 for each address of CHAMNHA #####
alert <- dplyr::select(all, alert_id, day, t2m, t2d, min_temp, max_temp, sd_temp, t2d_9, t2m_9 ,min_temp_9, max_temp_9, sd_temp_9, perc)
setDT(alert)
setkey(alert, alert_id, day)
gc()

### we estimate relative humidity
calculate_relative_humidity <- function(air_temp, dewpoint, m = 7.591386, 
                                        t_n = 240.7263) {
  # https://earthscience.stackexchange.com/questions/16570/how-to-calculate-relative-humidity-from-temperature-dew-point-and-pressure
  100 * 10 ^ (m * ((dewpoint / (dewpoint + t_n)) - (air_temp / (air_temp + t_n))))
}
alert[, rh := calculate_relative_humidity(t2m, t2d, m = 7.591386, t_n = 240.7263)]
alert[, rh_9 := calculate_relative_humidity(t2m_9, t2d_9, m = 7.591386, t_n = 240.7263)]

##### Create lag variables for temperature, up to lag 7
#     NOTE: I didn't use dplyr::mutate because the dataset is too large
## 28x28km
alert <- alert %>% group_by(alert_id) %>% mutate(temp_s1= dplyr::lag(t2m, n=1, default=NA),
                                                 temp_s2 = dplyr::lag(t2m, n=2, default=NA),
                                                 temp_s3 = dplyr::lag(t2m, n=3, default=NA),
                                                 temp_s4 = dplyr::lag(t2m, n=4, default=NA),
                                                 temp_s5 = dplyr::lag(t2m, n=5, default=NA),
                                                 temp_s6 = dplyr::lag(t2m, n=6, default=NA))

alert <- alert %>% group_by(alert_id) %>% mutate(rh_s1= dplyr::lag(rh, n=1, default=NA),
                                                 rh_s2 = dplyr::lag(rh, n=2, default=NA),
                                                 rh_s3 = dplyr::lag(rh, n=3, default=NA),
                                                 rh_s4 = dplyr::lag(rh, n=4, default=NA),
                                                 rh_s5 = dplyr::lag(rh, n=5, default=NA),
                                                 rh_s6 = dplyr::lag(rh, n=6, default=NA))

alert <- alert %>% group_by(alert_id) %>% mutate(temp_min_s1= dplyr::lag(min_temp, n=1, default=NA),
                                                 temp_min_s2 = dplyr::lag(min_temp, n=2, default=NA),
                                                 temp_min_s3 = dplyr::lag(min_temp, n=3, default=NA),
                                                 temp_min_s4 = dplyr::lag(min_temp, n=4, default=NA),
                                                 temp_min_s5 = dplyr::lag(min_temp, n=5, default=NA),
                                                 temp_min_s6 = dplyr::lag(min_temp, n=6, default=NA))

alert <- alert %>% group_by(alert_id) %>% mutate(temp_max_s1= dplyr::lag(max_temp, n=1, default=NA),
                                                 temp_max_s2 = dplyr::lag(max_temp, n=2, default=NA),
                                                 temp_max_s3 = dplyr::lag(max_temp, n=3, default=NA),
                                                 temp_max_s4 = dplyr::lag(max_temp, n=4, default=NA),
                                                 temp_max_s5 = dplyr::lag(max_temp, n=5, default=NA),
                                                 temp_max_s6 = dplyr::lag(max_temp, n=6, default=NA))

alert <- alert %>% group_by(alert_id) %>% mutate(temp_sd_s1= dplyr::lag(sd_temp, n=1, default=NA),
                                                 temp_sd_s2 = dplyr::lag(sd_temp, n=2, default=NA),
                                                 temp_sd_s3 = dplyr::lag(sd_temp, n=3, default=NA),
                                                 temp_sd_s4 = dplyr::lag(sd_temp, n=4, default=NA),
                                                 temp_sd_s5 = dplyr::lag(sd_temp, n=5, default=NA),
                                                 temp_sd_s6 = dplyr::lag(sd_temp, n=6, default=NA))

## 9x9km
alert <- alert %>% group_by(alert_id) %>% mutate(temp_s1_9= dplyr::lag(t2m_9, n=1, default=NA),
                                                 temp_s2_9 = dplyr::lag(t2m_9, n=2, default=NA),
                                                 temp_s3_9 = dplyr::lag(t2m_9, n=3, default=NA),
                                                 temp_s4_9 = dplyr::lag(t2m_9, n=4, default=NA),
                                                 temp_s5_9 = dplyr::lag(t2m_9, n=5, default=NA),
                                                 temp_s6_9 = dplyr::lag(t2m_9, n=6, default=NA))

alert <- alert %>% group_by(alert_id) %>% mutate(rh_s1_9 = dplyr::lag(rh_9, n=1, default=NA),
                                                 rh_s2_9 = dplyr::lag(rh_9, n=2, default=NA),
                                                 rh_s3_9 = dplyr::lag(rh_9, n=3, default=NA),
                                                 rh_s4_9 = dplyr::lag(rh_9, n=4, default=NA),
                                                 rh_s5_9 = dplyr::lag(rh_9, n=5, default=NA),
                                                 rh_s6_9 = dplyr::lag(rh_9, n=6, default=NA))

alert <- alert %>% group_by(alert_id) %>% mutate(temp_min_s1_9 = dplyr::lag(min_temp_9, n=1, default=NA),
                                                 temp_min_s2_9 = dplyr::lag(min_temp_9, n=2, default=NA),
                                                 temp_min_s3_9 = dplyr::lag(min_temp_9, n=3, default=NA),
                                                 temp_min_s4_9 = dplyr::lag(min_temp_9, n=4, default=NA),
                                                 temp_min_s5_9 = dplyr::lag(min_temp_9, n=5, default=NA),
                                                 temp_min_s6_9 = dplyr::lag(min_temp_9, n=6, default=NA))

alert <- alert %>% group_by(alert_id) %>% mutate(temp_max_s1_9 = dplyr::lag(max_temp_9, n=1, default=NA),
                                                 temp_max_s2_9 = dplyr::lag(max_temp_9, n=2, default=NA),
                                                 temp_max_s3_9 = dplyr::lag(max_temp_9, n=3, default=NA),
                                                 temp_max_s4_9 = dplyr::lag(max_temp_9, n=4, default=NA),
                                                 temp_max_s5_9 = dplyr::lag(max_temp_9, n=5, default=NA),
                                                 temp_max_s6_9 = dplyr::lag(max_temp_9, n=6, default=NA))

alert <- alert %>% group_by(alert_id) %>% mutate(temp_sd_s1_9 = dplyr::lag(sd_temp_9, n=1, default=NA),
                                                 temp_sd_s2_9 = dplyr::lag(sd_temp_9, n=2, default=NA),
                                                 temp_sd_s3_9 = dplyr::lag(sd_temp_9, n=3, default=NA),
                                                 temp_sd_s4_9 = dplyr::lag(sd_temp_9, n=4, default=NA),
                                                 temp_sd_s5_9 = dplyr::lag(sd_temp_9, n=5, default=NA),
                                                 temp_sd_s6_9 = dplyr::lag(sd_temp_9, n=6, default=NA))

alert <- alert %>% group_by(alert_id) %>% mutate(perc_s1= dplyr::lag(perc, n=1, default=NA),
                                                 perc_s2 = dplyr::lag(perc, n=2, default=NA),
                                                 perc_s3 = dplyr::lag(perc, n=3, default=NA),
                                                 perc_s4 = dplyr::lag(perc, n=4, default=NA),
                                                 perc_s5 = dplyr::lag(perc, n=5, default=NA),
                                                 perc_s6 = dplyr::lag(perc, n=6, default=NA))

setDT(alert)

# correct hospital data
table(alert$alert_id)

alert[, hosp := alert_id]
alert[, alert_id := NULL]
table(alert$hosp)
saveRDS(alert, "02_original/alert_temp_db.rds")

#### Work the the health data ####

##### I create a case-crossover dataset of the still birth and non-still birth
#al_hea <- read_dta("Appended - Feb 26 2024.dta")
#saveRDS(al_hea, "02_original/alert_health.rds")

al_hea <- readRDS("02_original/alert_health.rds")
setDT(al_hea)

# names
names(al_hea)

# hospital data
table(al_hea$hospn)
table(al_hea$hosp)

al_hea[, hospn := as.numeric(hospn)]
table(al_hea$hospn)
al_hea <- al_hea[!is.na(hospn)]

# country
table(al_hea$country)
al_hea[, country := as.factor(country)]
levels(al_hea$country) <- c("Benin", "Malawi", "Tanzania", "Uganda")

# q3age: maternal age
summary(al_hea$q3age)
hist(al_hea$q3age)

# q4grav: Number of pregnancy
summary(al_hea$q4grav)
hist(al_hea$q4grav)

# q9aga: gestational age
summary(al_hea$q9aga) # only 3000 missing data (2,1%)
hist(al_hea$q9aga)

# q12acom: multiple pregnancies
summary(al_hea$q12ccom)
table(al_hea$q12acom)
al_hea[,q12acom := as.factor(q12acom)]
levels(al_hea$q12acom) <- c("No", "Yes")

# q12ccom: preterm labour
summary(al_hea$q12ccom)
table(al_hea$q12ccom)
al_hea[,q12ccom := as.factor(q12ccom)]
levels(al_hea$q12ccom) <- c("No", "Yes")

# q23birth: date of birth
al_hea[, delivery_date := as.Date(as.POSIXct(q23birth, 'GMT'))]

# q24mode: mode of birth
table(al_hea$q24mode)
al_hea[q24mode == 1, deliv_mode := "Spontaneous"]
al_hea[q24mode == 2, deliv_mode := "Caesarean"]
al_hea[q24mode == 3 | q24mode == 4| q24mode == 5, deliv_mode := "Others"]
table(al_hea$deliv_mode)
al_hea$deliv_mode <- ordered(al_hea$deliv_mode, levels = c("Spontaneous", "Caesarean", "Others"))
summary(al_hea$deliv_mode)

# q28weight: birth weigght
summary(al_hea$q28weight)
hist(al_hea$q28weight)

# q29sex: sex
summary(al_hea$q29sex)
table(al_hea$q29sex)
al_hea[,q29sex := as.factor(q29sex)]
levels(al_hea$q29sex) <- c("Girl", "Boy")

# new variables: 
str(al_hea$q11hiv)
table(al_hea$q11hiv)
al_hea[,q11hiv := as.numeric(q11hiv)]
al_hea[q11hiv == 3, q11hiv := NA]
al_hea[q11hiv == 2, q11hiv := 0]
al_hea[,q11hiv := as.factor(q11hiv)]
levels(al_hea$q11hiv) <- c("No", "Yes")

str(al_hea$ht)
al_hea[,ht := as.factor(ht)]
levels(al_hea$ht) <- c("No", "Yes")

str(al_hea$prol_obst_lab)
al_hea[,prol_obst_lab := as.factor(prol_obst_lab)]
levels(al_hea$prol_obst_lab) <- c("No", "Yes")

str(al_hea$aph)
al_hea[,aph := as.factor(aph)]
levels(al_hea$aph) <- c("No", "Yes")

summary(al_hea$ht)
summary(al_hea$prol_obst_lab)
summary(al_hea$aph)

# stillbirth: stillbirth
table(al_hea$q27out)

# 1: No
# 2: Intrapartum stillbirth
# 3: Antepartum stillbirth

l <- lapply(al_hea, attr, "label")
l
#NEW
# create subset:
# intrastill birth
al_hea[q27out == 2, intra_still := 1]
al_hea[q27out == 1, intra_still := 0]
al_hea[q27out == 3, intra_still := 0]

# antepartum birth
al_hea[q27out == 3, ante_still := 1]
al_hea[q27out == 2, ante_still := 0]
al_hea[q27out == 1, ante_still := 0]

# all strill birth
al_hea[q27out == 3, all_still := 1]
al_hea[q27out == 2, all_still := 1]
al_hea[q27out == 1, all_still := 0]

# perinatal deaths (stillbirths + neonatal death)
al_hea[all_still == 1 | neo_death == 1, per_mort :=1]
al_hea[all_still == 0 & neo_death == 0, per_mort := 0]


table(al_hea$q27out)
table(al_hea$intra_still)
table(al_hea$ante_still)
table(al_hea$all_still)
table(al_hea$per_mort)
al_hea[is.na(per_mort) & !is.na(all_still), .(q27out, all_still, intra_still,ante_still,neo_death, per_mort)]

# neo_death: neonatal death
table(al_hea$neo_death)

#### convert
al_hea[,all_still := as.factor(all_still)]
al_hea[,intra_still := as.factor(intra_still)]
al_hea[,ante_still := as.factor(ante_still)]
al_hea[,per_mort := as.factor(per_mort)]
al_hea[,neo_death := as.factor(neo_death)]

levels(al_hea$all_still) <- c("No", "Yes")
levels(al_hea$intra_still) <- c("No", "Yes")
levels(al_hea$ante_still) <- c("No", "Yes")
levels(al_hea$neo_death) <- c("No", "Yes")
levels(al_hea$per_mort) <- c("No", "Yes")

#### remove children who have both missing on neo_death and others
al_hea[is.na(q27out) & is.na(neo_death), rem := TRUE]
al_hea[rem == TRUE]

al_hea <- al_hea[is.na(rem)]
al_hea[,rem := NULL]

# check if there are some stillbirth and neonataldeath
al_hea[q27out >1 & neo_death == 1, rem := TRUE]
al_hea[rem == TRUE] # no overlaå
al_hea[,rem := NULL]

# if there was no still birth but information of neonatal death consider them as 1 and neonatal deaths
al_hea[per_mort == 1 & is.na(all_still), all_still := 0]
al_hea[per_mort == 1 & is.na(intra_still), intra_still := 0]
al_hea[per_mort == 1 & is.na(ante_still), ante_still := 0]

al_hea[per_mort == 1 & is.na(q27out), .(q27out, all_still, intra_still, ante_still,per_mort)]

# subset the data with for creating the subset variables
pt <- al_hea[,.(alert_id,delivery_date )]

# create alert dataset with relevant variables
alert_db <- al_hea[,.(alert_id, hospn, country, q3age, q4grav, q9aga, q12acom, q12ccom, deliv_mode,
                      q23birth, q28weight, q29sex, q27out, q11hiv, prol_obst_lab, ht, aph, intra_still, ante_still, all_still, per_mort, neo_death,referred )]

# save clean alert_db
saveRDS(alert_db, "02_original/alert_health_db.rds")

##### CREATE CC dataset #####

library(rSPARCS)
library(dplyr)
pt <- as.data.frame(pt)

# match each case day with same weekdays in the same month  
cc1 <- CXover.data(data=pt, date="delivery_date", ID="alert_id")
setDT(cc1)

## delete the first column which is the newly created ID for participants
cc1[,1] <- NULL
setnames(cc1, "Date","day")
setDT(cc1)

# for know we eliminate date before 2 years and any missing values
alert <- alert[!is.na(temp_s6)]  # remove missing

# merge with health data
cc1 <- merge(cc1, alert_db, by = "alert_id", all.x = TRUE)

##### Finally I merge the CC file with the (lagged) temperatures file, and save as CC dataset
setnames(cc1, "hospn", "hosp") # change the names,,
cc2 <- merge(cc1, alert, by = c("hosp", "day"), all.x = TRUE)
cc2 <- cc2[order(alert_id, day)]

#change the names
setnames(cc2, c("t2m","rh","min_temp","max_temp","sd_temp","t2m_9","t2d_9","min_temp_9","max_temp_9","sd_temp_9"),
         c("temp_s0","rh_s0", "temp_min_s0","temp_max_s0","temp_sd_s0","temp_s0_9","rh_s0_9", "temp_min_s0_9","temp_max_s0_9","temp_sd_s0_9"))

# hospital is a factor
cc2[, hosp := as.factor(hosp)]

############ check full data or not

# for now we don´t have data for full mionth of may! so exclude if date of birth >01
al_hea <- readRDS("02_original/alert_health.rds")
length(unique(al_hea$alert_id)) # N = 138213

# remove that does not have data on perinatal death
cc2[is.na(intra_still)]
rem  <- cc2[is.na(per_mort)]
cc2  <- cc2[!is.na(per_mort)]

# check numbers
length(unique(rem$alert_id)) # 51 + 30 that previously were excluded = 81
length(unique(cc2$alert_id)) # 138160

# remove exposure data
rem  <- cc2[is.na(temp_s0)] 
cc2  <- cc2[!is.na(temp_s0)]

# check numbers
length(unique(rem$alert_id)) # 145
length(unique(cc2$alert_id)) # 138015
 
saveRDS(cc2, "03_data/alert_cc.rds")

##### create longitudinal database ####

## open temperature db
temp <- readRDS("02_original/alert_temp_db.rds")

## health dataset
alert_db <- readRDS("02_original/alert_health_db.rds")

## We are going to estimate from the date of birth and gestational age a weekly date (going back).
## once we have the date we can link the temperature exposures
library(tidyverse)

db <- alert_db[,.(alert_id, hospn, q23birth, q9aga)]
colnames(db) <- c("alert_id","hosp" ,"date_of_birth", "gestational_age")
db[, date_of_birth := as.Date(date_of_birth, format= "%Y-%m-%d")]
db <- db[!is.na(gestational_age)]

# Create a function to generate a sequence of gestational weeks in a retrospective way
generate_gestational_weeks <- function(gestational_age) {
  if (!is.finite(gestational_age)) gestational_age <- NA
  seq(1, gestational_age, by = 1)
}

# Apply the function to create a long-format dataset
long_data <- db %>%
  rowwise() %>%
  mutate(gestational_week = list(generate_gestational_weeks(gestational_age))) %>%
  unnest(cols = gestational_week) %>%
  mutate(date = date_of_birth - (gestational_age - gestational_week) * 7)

setDT(long_data)
hist(long_data$gestational_week)
table(long_data$gestational_week)

# extract season of conception
ges_wk <- long_data[ gestational_week == 1, .(alert_id, day )]
ges_wk[, con_m := as.numeric(format(day, "%m"))]
ges_wk[, con_y := as.numeric(format(day, "%Y"))]

## looks perfect. Now we can link the exposure database
setnames(long_data, "date", "day")
long <- merge(long_data, temp, by = c("hosp", "day"), all.x = TRUE) # no we have the daily exposure for each day of the week

# change names
setnames(long, c("t2m","rh","min_temp","max_temp","sd_temp","t2m_9","th_9","min_temp_9","max_temp_9","sd_temp_9"),
         c("temp_s0","rh_s0", "temp_min_s0","temp_max_s0","temp_sd_s0","temp_s0_9","th_s0_9", "temp_min_s0_9","temp_max_s0_9","temp_sd_s0_9"))

# estimate the weekly levels 28x28
long[, temp_wk := (temp_s0 + temp_s1 + temp_s2 + temp_s3 + temp_s4 + temp_s5 + temp_s6)/7]
long[, rh_wk := (rh_s0 + rh_s1 + rh_s2 + rh_s3 + rh_s4 + rh_s5 + rh_s6)/7]
long[, min_wk := (temp_min_s0 + temp_min_s1 + temp_min_s2 + temp_min_s3 + temp_min_s4 + temp_min_s5 + temp_min_s6)/7]
long[, max_wk := (temp_max_s0 + temp_max_s1 + temp_max_s2 + temp_max_s3 + temp_max_s4 + temp_max_s5 + temp_max_s6)/7]
long[, sd_wk := (temp_sd_s0 + temp_sd_s1 + temp_sd_s2 + temp_sd_s3 + temp_sd_s4 + temp_sd_s5 + temp_sd_s6)/7]

# estimate the weekly levels 9x9km
long[, temp_wk_9 := (temp_s0_9 + temp_s1_9 + temp_s2_9 + temp_s3_9 + temp_s4_9 + temp_s5_9 + temp_s6_9)/7]
long[, rh_wk_9 := (rh_s0_9 + rh_s1_9 + rh_s2_9 + rh_s3_9 + rh_s4_9 + rh_s5_9 + rh_s6_9)/7]
long[, min_wk_9 := (temp_min_s0_9 + temp_min_s1_9 + temp_min_s2_9 + temp_min_s3_9 + temp_min_s4_9 + temp_min_s5_9 + temp_min_s6_9)/7]
long[, max_wk_9 := (temp_max_s0_9 + temp_max_s1_9 + temp_max_s2_9 + temp_max_s3_9 + temp_max_s4_9 + temp_max_s5_9 + temp_max_s6_9)/7]
long[, sd_wk_9 := (temp_sd_s0_9 + temp_sd_s1_9 + temp_sd_s2_9 + temp_sd_s3_9 + temp_sd_s4_9 + temp_sd_s5_9 + temp_sd_s6_9)/7]

# identify the periods of week
long[gestational_week >= 0 & gestational_week <13, trimester := "first"]
long[gestational_week >= 13 & gestational_week <28, trimester := "second"]
long[gestational_week >= 28, trimester := "third"]
p <-long[,.(gestational_week, trimester)]
p #looks correct

# estimate gestational week using inverse format. 
long[, inv_gest := gestational_age - gestational_week]
p <- long[, .(inv_gest, gestational_age, gestational_week)]
p # week 0 is the week of birth, week 40 is the closest week of conception

# estimate the mean levels of temperature during pregnancy
temp_preg <- long %>%
  group_by(alert_id) %>%
  summarize(temp_preg = mean(temp_wk, na.rm = TRUE),
            rh_preg = mean(rh_wk, na.rm = TRUE),
            min_preg = mean(min_wk, na.rm = TRUE),
            max_preg = mean(max_wk, na.rm = TRUE),
            sd_preg = mean(sd_wk, na.rm = TRUE),
            temp_preg_9 = mean(temp_wk_9, na.rm = TRUE),
            rh_preg_9 = mean(rh_wk_9, na.rm = TRUE),
            min_preg_9 = mean(min_wk_9, na.rm = TRUE),
            max_preg_9 = mean(max_wk_9, na.rm = TRUE),
            sd_preg_9 = mean(sd_wk_9, na.rm = TRUE))


# estimate the levels of temperature for each trimester
temp_trim <- long %>%
  group_by(alert_id, trimester) %>%
  summarize(temp_wk = mean(temp_wk, na.rm = TRUE),
            rh_wk = mean(rh_wk, na.rm = TRUE),
            min_wk = mean(min_wk, na.rm = TRUE),
            max_wk = mean(max_wk, na.rm = TRUE),
            sd_wk = mean(sd_wk, na.rm = TRUE),
            temp_wk_9 = mean(temp_wk_9, na.rm = TRUE),
            rh_wk_9 = mean(rh_wk_9, na.rm = TRUE),
            min_wk_9 = mean(min_wk_9, na.rm = TRUE),
            max_wk_9 = mean(max_wk_9, na.rm = TRUE),
            sd_wk_9 = mean(sd_wk_9, na.rm = TRUE))
setDT(temp_trim)

# export to wide format
trim <- dcast(melt(temp_trim, id.vars=c("alert_id", "trimester")), alert_id~variable+trimester)

# merge with preg and trimester exposure
alert_long <- merge(temp_preg, trim, by = "alert_id", all = TRUE) 
# this data includes for each temperature metric the mean exposure during pregnancy and by trimester

## we are now going to select the weekly variables from the gestational age week (for now only weekly exposures)
week_exp <- long[,.(alert_id, inv_gest, temp_wk , temp_wk_9)]
colnames(week_exp) <- c("alert_id", "inv_gest", "temp_week" , "temp_week_9")
week_exp <- dcast(melt(week_exp, id.vars=c("alert_id", "inv_gest")), alert_id~variable+inv_gest)

# merge with original file
alert_long <- merge(alert_long, week_exp, by = "alert_id", all = TRUE)

names(alert_long)

# temp_preg := mean temperature during pregnancy
# temp_wk_third := trimester specific temperatures
# temp_week := weekly mean temperatures for each gestational week.

## add the health variables but firt the weeks of conceptions
alert_db

# merge weeks of conception
ges_wk <- ges_wk[,.(alert_id, con_m, con_y)]
alert <- merge(alert_db, ges_wk, by = "alert_id", all = TRUE)

# merge exposures
alert <- merge(alert, alert_long, by = "alert_id", all = TRUE)

## check duplicates
#by ID
alert[, dup := 1:.N, by = "alert_id"]
table(alert$dup) # no duplicates
alert[, dup := NULL]

### people without outcome data
summary(alert$all_still)
summary(alert$ante_still)
summary(alert$intra_still)
summary(alert$per_mort) # few missing. Remove missing otucome data

# remove missing outcome data
db_alert <- alert[!is.na(all_still)| !is.na(per_mort)] 
rem  <- alert[is.na(all_still)| is.na(per_mort)] 
length(unique(rem$alert_id)) # 71 exclusions
length(unique(db_alert$alert_id)) # 138214 unique data

# remove missing exposure data
rem  <- db_alert[is.na(temp_preg)] 
length(unique(rem$alert_id)) # 3005 missing

# 2900 are because no gestational age
rem  <- rem[!is.na(q9aga)] # the other 100 don´t have data on birthday

# remove
db_alert <- db_alert[!is.na(temp_preg)] 
length(unique(db_alert$alert_id)) # 135209 unique data


## save
saveRDS(db_alert, "03_data/alert_long.rds")





