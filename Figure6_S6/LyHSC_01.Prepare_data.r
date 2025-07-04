# '''
# This is to analyse HSC1 and HSC2.
# '''
rm(list=ls())
set.seed(1234)

print('---------------------------------')
print(Sys.Date())
print('---------------------------------')


source("/home/mayao/script/Function.r")



### ------------------------------------------------
###                   Function Zone
### ------------------------------------------------
# findMarkers2Groups: find markers and padj cut and decreasing order by log2FC.
findMarkers2Groups <- function(obj, ident1, ident2){
    marker <- FindMarkers(hsc, ident.1 = ident1, ident.2 = ident2)  
    marker <- marker[marker$p_val_adj < 0.05,]
    marker <- marker[order(marker$avg_log2FC, decreasing = T),]
}


### ------------------------------------------------
###                   Running Zone
### ------------------------------------------------

### --------------------
### 01. Set up env
### --------------------
suppressMessages(library(dplyr))
library(ggplot2) %>% suppressMessages()
library(Seurat) %>% suppressMessages()
library(plyr) %>% suppressMessages()
library(cowplot) %>% suppressMessages()
library(sctransform) %>% suppressMessages()
library(stringr) %>% suppressMessages()
# library(harmony) %>% suppressMessages()
library(scales) %>% suppressMessages()
### set path
setwd("/home/mayao/LABdata/MY/hsc1_hsc2/mpp")


### --------------------
### 02. Prepare data
### --------------------
## read in data
meta <- openxlsx::read.xlsx("meta/I05.MPP.xlsx", sheet = 1, rowNames = FALSE)
rownames(meta) <- meta$Cell

## remove non-coding and duplicated genes
origData <- read.table("count/MPP.rawCount.txt", header = T) # 53379 291
cleanData <- keepCodingGenes(origData) %>% removeDuplicatedGenes()  # 21940 289
## only coding genes and no duplication
# column names: lib1_sc01 lib1_sc02 ... lib8_sc96 ENSEMBL
cleanMeta <- meta[colnames(cleanData[1:ncol(cleanData)-1]) ,] 
cleanMeta <- cleanMeta[!is.na(cleanMeta$Group),]
print(cleanMeta[1:30,])
cleanData <- cleanData[, c(rownames(cleanMeta),"ENSEMBL")]
print(dim(cleanData)) # 21940   277
print(dim(cleanMeta)) # 276   7

## save cleanData and cleanMeta
save(cleanData, cleanMeta, file = "Rdata/01.mpp.cleanData.cleanMeta.rdata" )


### ----------------------
### 03. Create Seurat obj
### ----------------------
load("Rdata/01.mpp.cleanData.cleanMeta.rdata")
dim(cleanData)
mpp <- CreateSeuratObject(counts = cleanData[,1:ncol(cleanData)-1], 
                                  meta.data = cleanMeta, 
                                  assay = "RNA", 
                                  min.cells = 3, 
                                  min.features = 200)
mpp <- PercentageFeatureSet(mpp, 
                           pattern = "^mt-", 
                           col.name = "pct_counts_Mito")
dim(mpp) # 11686   256
table(mpp$Group)
# MPP5 MPP6 MPP7 
#   86   84   86 
is.na(mpp$CD135_median) %>% table()
# FALSE  TRUE 
#   250     6

save(mpp,file = 'Rdata/temp01.mpp.seuratObj.rdata')
