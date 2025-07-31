#!/usr/bin/env Rscript
################################################################################
# AP-DNAme-LF - Test for Gene Set Enrichment (CpGs + AMRs)
################################################################################
# Author: Lucile
# Date of creation:20/03/2023
# Note: 
# Run Rscript 3-HDMAX2.R before
################################################################################
# Paths and parameters
#------------------------------------------------------------------------------#

# All DMRs: 
resDir <- "hdmax2"

saveDir <- "GSEA"

probevarFile <- "Data/SEP_ProbeVariation_PRadj.rds"

all.cpg <- rownames(readRDS(file = "/Data/DNAmethylplacenta_adjfunnorm_preprocessed.rds"))
length(all.cpg)
# 755364
# Parameters (GSEA)

arraytype <- "EPIC"
ncores <- 2

################################################################################
# R packages
#------------------------------------------------------------------------------#

library(magrittr)
library(ggplot2)
library(missMethyl)
library(GenomicRanges)

suppressPackageStartupMessages( library( IlluminaHumanMethylationEPICanno.ilm10b4.hg19 ) ) 
annotation_object <- IlluminaHumanMethylationEPICanno.ilm10b4.hg19::IlluminaHumanMethylationEPICanno.ilm10b4.hg19
FullAnnot <- minfi::getAnnotation(annotation_object)

################################################################################
# R function(s)
#------------------------------------------------------------------------------#

source("Rscripts/Rfunctions.R")

# Credit to Alexandra Binder

createCollectionFromDB <- function(myDBfile){
  
  library(org.Hs.eg.db)
  
  # External updated database
  dat <- scan(myDBfile, what="", sep="\n")
  ldat <- strsplit(dat, "\t")
  names <- sapply(ldat,function(x)x[1])
  ldat <- lapply(ldat, function(x) x[-1])
  ldat <- lapply(ldat, function(x) x[x != ""])
  names(ldat)<-names
  #class(ldat)
  
  ### Converting Gene Symbols to Entrez ID
  ldat_entrez<-lapply(ldat,function(x){
    
    suppressMessages(AnnotationDbi::select(org.Hs.eg.db, 
                                           keys = x,
                                           columns = c("ENTREZID", "SYMBOL"),
                                           keytype = "SYMBOL")$ENTREZID)
    
  })
  
  return( ldat_entrez )
}

################################################################################
# Probe variation
#------------------------------------------------------------------------------#

probeVar <- readRDS("Data/SEP_ProbeVariation_PRadj.rds")

# Background (all *tested* probes)

inFile <- "Hits_AMRs.rds"
dmrs <- readRDS(inFile)

inFile.cpgs <- "Hits_CpG.xlsx"
dmps <- xlsx::read.xlsx(inFile.cpgs, sheetIndex = 1)

if( !dir.exists(saveDir) ) dir.create( saveDir )

exposure.names <- unique(dmrs$Exposure)
outcome.names <- unique(dmrs$Outcome)

################################################################################
# RUN
#------------------------------------------------------------------------------#

ldat_entrez <- createCollectionFromDB(myDBfile = "Databases/KEGG_2021_Human")

for(o in seq_along( outcome.names )){
  
  GSEA <- data.frame()
  
  for(f in seq_along(exposure.names)){
    
    cat(paste0(exposure.names[f], " - ", outcome.names[o], "\n"))
    
    res <- dmrs %>% dplyr::filter(Outcome == outcome.names[o] & Exposure == exposure.names[f])
    res <- res %>% dplyr::filter(!is.na(geneName) & geneName != "")
    
    res2 <- dmps %>% dplyr::filter(Outcome.org == outcome.names[o] & Exposure.org == exposure.names[f]) 
    res2 <- res2 %>% dplyr::filter(!is.na(geneName) & geneName != "")
    
    dim(res) %>% print()

  if( nrow(res)>0 ){
  # regions: GRanges object of DMR coordinates to test for GO term enrichment
  regions <- GRanges(seqnames = res$DMR.chr, ranges = IRanges(start = res$DMR.start,end = res$DMR.end))
  sites <- GRanges(seqnames = res2$chr, ranges = IRanges(start = res2$pos,end = res2$pos))
  
  regions <- union(regions, sites)
  
  # KEGG pathway enrichment
  
  # Enrichment 
  df_tmp <- try(missMethyl::gsaregion(regions = regions, 
                                     all.cpg = all.cpg,
                                     collection = ldat_entrez,
                                     prior.prob = TRUE,
                                     anno = FullAnnot,
                                     genomic.features = c("ALL"),
                                     array.type = arraytype, 
                                     sig.genes = T))
  df_tmp <- df_tmp %>% dplyr::filter(N>20)
  df_tmp$FDR2 <- p.adjust(p = df_tmp$P.DE, method = "BH")
  df_tmp <- cbind.data.frame(Name = rownames(df_tmp), df_tmp)
  
  if( class(df_tmp) != "try-error" ){
  if( nrow(df_tmp) > 0 ) GSEA <- rbind.data.frame(GSEA,
                                                  data.frame(Exposure = exposure.names[f],
                                                             Outcome = outcome.names[o],
                                                             geneSet = "all",
                                                             location = "ALL",
                                                             df_tmp))
  
  

  
    }# End try-error
    }# End if results
    outFile <- paste0(saveDir, "/KEGG_regions.dmrs_",outcome.names[o],".rds")
    
    saveRDS(GSEA, file = outFile)
    
    GSEA %>% dplyr::filter(FDR2<0.2) %>% print()
    #GSEA %>% dplyr::arrange(P.DE) %>% head(10) %>% print()
  }# End loop on exposures
}# End loop on outcomes

################################################################################
# Reactome
#------------------------------------------------------------------------------#

ldat_entrez <- createCollectionFromDB(myDBfile = "/Databases/Reactome_2022")

for(o in seq_along( outcome.names )){
  
  GSEA <- data.frame()
  
  for(f in seq_along(exposure.names)){
    
    cat(paste0(exposure.names[f], " - ", outcome.names[o], "\n"))
    
    res <- dmrs %>% dplyr::filter(Outcome == outcome.names[o] & Exposure == exposure.names[f])
    res <- res %>% dplyr::filter(!is.na(geneName) & geneName != "")
    
    res2 <- dmps %>% dplyr::filter(Outcome.org == outcome.names[o] & Exposure.org == exposure.names[f]) 
    res2 <- res2 %>% dplyr::filter(!is.na(geneName) & geneName != "")
    
    dim(res) %>% print()
    
    if( nrow(res)>0 ){
      # regions: GRanges object of DMR coordinates to test for GO term enrichment
      regions <- GRanges(seqnames = res$DMR.chr, ranges = IRanges(start = res$DMR.start,end = res$DMR.end))
      sites <- GRanges(seqnames = res2$chr, ranges = IRanges(start = res2$pos,end = res2$pos))
      
      regions <- union(regions, sites)
      
      # KEGG pathway enrichment
      
      # Enrichment 
      df_tmp <- try(missMethyl::gsaregion(regions = regions, 
                                          all.cpg = all.cpg,
                                          collection = ldat_entrez,
                                          prior.prob = TRUE,
                                          anno = FullAnnot,
                                          genomic.features = c("ALL"),
                                          array.type = arraytype, 
                                          sig.genes = T))
      df_tmp <- df_tmp %>% dplyr::filter(N>20)
      df_tmp$FDR2 <- p.adjust(p = df_tmp$P.DE, method = "BH")
      df_tmp <- cbind.data.frame(Name = rownames(df_tmp), df_tmp)
      
      if( class(df_tmp) != "try-error" ){
        if( nrow(df_tmp) > 0 ) GSEA <- rbind.data.frame(GSEA,
                                                        data.frame(Exposure = exposure.names[f],
                                                                   Outcome = outcome.names[o],
                                                                   geneSet = "all",
                                                                   location = "ALL",
                                                                   df_tmp))
        
        
        
        
      }# End try-error
    }# End if results
    outFile <- paste0(saveDir, "/Reactome_",outcome.names[o],".rds")
    
    saveRDS(GSEA, file = outFile)
    
    GSEA %>% dplyr::filter(FDR2<0.2) %>% print()
    #GSEA %>% dplyr::arrange(P.DE) %>% head(10) %>% print()
  }# End loop on exposures
}# End loop on outcomes

################################################################################
# Go term enrichment
#------------------------------------------------------------------------------#

ldat_entrez <- createCollectionFromDB(myDBfile = "/Databases/GO_Biological_Process_2023")

for(o in seq_along( outcome.names )){
  
  GSEA <- data.frame()
  
  for(f in seq_along(exposure.names)){
    
    cat(paste0(exposure.names[f], " - ", outcome.names[o], "\n"))
    
    res <- dmrs %>% dplyr::filter(Outcome == outcome.names[o] & Exposure == exposure.names[f])
    res <- res %>% dplyr::filter(!is.na(geneName) & geneName != "")
    
    res2 <- dmps %>% dplyr::filter(Outcome.org == outcome.names[o] & Exposure.org == exposure.names[f]) 
    res2 <- res2 %>% dplyr::filter(!is.na(geneName) & geneName != "")
    
    dim(res) %>% print()
    
    if( nrow(res)>0 ){
      # regions: GRanges object of DMR coordinates to test for GO term enrichment
      regions <- GRanges(seqnames = res$DMR.chr, ranges = IRanges(start = res$DMR.start,end = res$DMR.end))
      sites <- GRanges(seqnames = res2$chr, ranges = IRanges(start = res2$pos,end = res2$pos))
      
      regions <- union(regions, sites)
      
      # KEGG pathway enrichment
      
      # Enrichment 
      df_tmp <- try(missMethyl::gsaregion(regions = regions, 
                                          all.cpg = all.cpg,
                                          collection = ldat_entrez,
                                          prior.prob = TRUE,
                                          anno = FullAnnot,
                                          genomic.features = c("ALL"),
                                          array.type = arraytype, 
                                          sig.genes = T))
      df_tmp <- df_tmp %>% dplyr::filter(N>20)
      df_tmp$FDR2 <- p.adjust(p = df_tmp$P.DE, method = "BH")
      df_tmp <- cbind.data.frame(Name = rownames(df_tmp), df_tmp)
      
      if( class(df_tmp) != "try-error" ){
        if( nrow(df_tmp) > 0 ) GSEA <- rbind.data.frame(GSEA,
                                                        data.frame(Exposure = exposure.names[f],
                                                                   Outcome = outcome.names[o],
                                                                   geneSet = "all",
                                                                   location = "ALL",
                                                                   df_tmp))
        
        
        
        
      }# End try-error
    }# End if results
    outFile <- paste0(saveDir, "/GOBP_",outcome.names[o],".rds")
    
    saveRDS(GSEA, file = outFile)
    
    GSEA %>% dplyr::filter(FDR2<0.2) %>% print()
    #GSEA %>% dplyr::arrange(P.DE) %>% head(10) %>% print()
  }# End loop on exposures
}# End loop on outcomes

################################################################################
