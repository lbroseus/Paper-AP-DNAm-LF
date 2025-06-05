#!/usr/bin/env Rscript
rm(list = ls()); gc()
################################################################################
# SEPAGES-DNAme - Preprocessing methylation data from norm. beta values
################################################################################
# Author: Lucile
# Date: 03/09/2024
# Update: 04/11/2024
################################################################################
# Steps:
# 1. Normalisation: InterpolatedXY adjusted funnorm (already performed)
# 2. Removal of CpGs close to known SNPs
# 3. Removal of cross-reactive CpGs 
# 4. Offset 0 and 1 beta values
# 5. Transformation to M-values
# --> Save meth data set for QC and technical var estimation 
# 6. Imputation of outlying probes (1% winsorizing)
# (cf: https://www.epicenteredresearch.com/pace/birthsize/)
# --> Save meth datasets fit for EWAS analyses
#------------------------------------------------------------------------------#
# Notes:
################################################################################
# Install R packages
#------------------------------------------------------------------------------#

# Available on the CRAN
install.packages(c("magrittr", "DescTools"))

# Available via Bioconductor
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("minfi")
BiocManager::install("DMRcate")
# necessary for minfi
BiocManager::install("bumphunter") 
BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19") # necessary for DMRcate
BiocManager::install("minfiData") # necessary for maxprobes

# Via Github
install.packages("remotes")
remotes::install_github("markgene/maxprobes")

################################################################################
# Paths and parameters
#------------------------------------------------------------------------------#

# Normalized data
methFile <- "DNAmethylplacenta_adjfunnorm.rds"

array_type <- "EPIC"

saveDir <- "Data/"

# Final datasets, ready for EWAS on autosomes 
outFile <- "DNAmethylplacenta_adjfunnorm_preprocessed.rds"

################################################################################
# R packages
#------------------------------------------------------------------------------#

library(magrittr)
library(DMRcate)
library(maxprobes)
#library(minfi)
#library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

################################################################################
# R functions
#------------------------------------------------------------------------------#

correctBorderlineValues <- function(probe){
  
  # 0 -> offset for null intensities
  if(min(probe, na.rm = T)==0) probe[probe==0] <- min(probe[probe>0])/10
  if(max(probe, na.rm = T)==1) probe[probe==1] <- max(probe[probe<1])+(1-max(probe[probe<1]))/10
  
  return( probe )
}

windsorizeOutliers <- function(probe, percent){
  
  if(percent<0 | percent>0.15) stop("percent must be positive and below 15%\n")
  
  probe <- DescTools::Winsorize(probe, val = quantile(probe, probs = c(percent/2, 1-percent/2)))
  
  return( probe )
}

################################################################################
# 1. Get normalized beta values (custom) and design
#------------------------------------------------------------------------------#

meth <- readRDS(methFile)

dim(meth)
# 824979    395

if(nrow(meth) < ncol(meth)) meth <- t( meth )

# check these are beta values
stopifnot(min(meth, na.rm = T)>=0 & max(meth,na.rm = T)<=1)

################################################################################
# 2. Removal of CpGs close to SNPs (dist<2 bp) and CHR X and Y
#------------------------------------------------------------------------------#

meth <- DMRcate::rmSNPandCH(meth, 
                            dist = 2, 
                            mafcut = 0.05,
                            rmcrosshyb = FALSE, 
                            rmXY = TRUE)

################################################################################
# 3. Removal of known cross-reactive probes
#------------------------------------------------------------------------------#

xloci <- maxprobes::xreactive_probes(array_type = array_type)
length(xloci)
# 43256 known cross-reactive probes in EPIC arrays

x <- which(rownames(meth) %in% xloci)
length(x)
# 40396 known cross-reactive probes still in the dataset

meth <- meth[-x,]

################################################################################
# 4. Imput offset value for intensities equal to 0 or 1 
# (required for M-value transformation)
#------------------------------------------------------------------------------#

meth <- apply(meth, 1, function(x) correctBorderlineValues(x))

meth <- t(meth)

################################################################################
# 5. Transformation into M-values
# https://rdrr.io/github/xuz1/ENmix/man/B2M.html
#------------------------------------------------------------------------------#

meth <- minfi::logit2(meth)

if(nrow(meth)<ncol(meth)) meth <- t(meth)

################################################################################
# 6. Winsorizing of outlying probes (1% of extreme values; 0.05/0.95)
#------------------------------------------------------------------------------#

meth <- apply(meth, 1, function(x) windsorizeOutliers(x, percent = 0.01))

meth <- t( meth )

################################################################################
# Save
#------------------------------------------------------------------------------#

saveRDS(meth, file = outFile)

################################################################################
