rm(list = ls())
################################################################################
# AirPollution-DNAm-LF - [2m] impute missing values (clinical covariates) 
################################################################################
# Author: Lucile
# Date: 19/01/2023 
# Notes: 
# -
################################################################################
# Covariates:
# ChildSex
# MaternalAge
# Parity
# MaternalBMI
# ParentalEducation
# MaternalSmoking
# ChildPassiveSmoking (12m)
# EDI
################################################################################
# R packages
#------------------------------------------------------------------------------#

library(magrittr)

################################################################################
# 
#------------------------------------------------------------------------------#

metaFile <- "SEP_metaData_20240617.rds"

expoFile <- "SEP_AirPollutionData.rds"

metaFile.imp <- "SEP_metaData_imp_2m.rds"

################################################################################
# Set seed for imputation
#------------------------------------------------------------------------------#

set.seed(38)

################################################################################
# Load metadata
#------------------------------------------------------------------------------#

metaData <- readRDS(metaFile)
dim(metaData)
# 484  51

metaData <- merge(metaData, readRDS(expoFile), all.x = T)
dim(metaData)
# 484  87

################################################################################
# Load metadata
#------------------------------------------------------------------------------#

outcomes <- c("LCI", "FRC", "TV", "RR", "MinVent", "MVperWe", "tPTEF_tE")
exposures <- c("NO2_p", "PM25_p", "PM10_p")

continuous_covariates <- c("MaternalAge", 
                           "MaternalSmoking_avgr",
                           "MaternalBMI", 
                           "EDI",
                           "ChildLength.2m",
                           "ChildWeight.2m", 
                           "Breastfeeding",
                           "ChildAge.2m")

categorical_covariates <- c("ChildSex", 
                            "MaternalSmoking",
                            "Ethnicity",
                            "Parity",
                            "ParentalEducation",
                            "MaternalEducation",
                            "ParentalRhinitis",
                            "DeliveryMode",
                            "ChildPassiveSmoking.12m",
                            "SeasonExam.2m", 
                            "SeasonBirth", 
                            "SeasonConception")

################################################################################
# Impute missing data 
#------------------------------------------------------------------------------#

metaData.imp <- metaData[,c("id", outcomes,exposures,continuous_covariates,categorical_covariates)]
dim(metaData.imp)
#395  20

bounds <- matrix(data=c(16,0,48,   #Breastfeeding must be >=0 and <=48 weeks
                        17,30,120, # Bound child age at the exam (default)
                        12,0,50  #MaternalSmoking must be >=0 and <50
                        ), 
                 ncol = 3, byrow = T)
amelia_fit <- Amelia::amelia(metaData.imp[,-1], 
                             m = 1, 
                             parallel = "multicore", 
                             noms = categorical_covariates,
                             bounds = bounds)

amelia_fit <- amelia_fit$imputations[[1]]
amelia_fit <- amelia_fit[,c(continuous_covariates,categorical_covariates)]

metaData.imp[,c(continuous_covariates,categorical_covariates)] <- amelia_fit

################################################################################
# Save
#------------------------------------------------------------------------------#

saveRDS(metaData.imp, file = metaFile.imp)

################################################################################
