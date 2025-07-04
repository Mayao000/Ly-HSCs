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
setwd("/home/mayao/LABdata/MY/hsc1_hsc2/new")


### --------------------
### 02. Prepare data
### --------------------
## read in data
meta <- read.table("meta/HSC1_HSC2_my.meta.txt", header = T)
rownames(meta) <- meta$Cell
# meta$Library <- sapply(meta$Library, function(x){
#     gsub('lib','Library',x)
# })
# meta$Cell <- sapply(meta$Cell, function(x){
#     gsub('lib','Library',x)
# })

# rownames(meta) <- meta$Cell
# meta$Group <- factor(meta$Group, levels = c('HSC1', 'HSC2'))
# write.table(meta, file='meta/HSC1_HSC2_my.meta.txt', sep = '\t', row.names = F, quote = F)

## remove non-coding and duplicated genes
origData <- read.table("matrix/HSC1_HSC2_my.rawCount.txt", header = T)
cleanData <- keepCodingGenes(origData) %>% removeDuplicatedGenes()  ## only coding genes and no duplication
                                                                    # column names: lib1_sc01 lib1_sc02 ... lib8_sc96 ENSEMBL
cleanMeta <- meta[colnames(cleanData[1:ncol(cleanData)-1]) ,] 
cleanMeta <- cleanMeta[!is.na(cleanMeta$Group),]
print(cleanMeta[1:30,])
cleanData <- cleanData[, c(rownames(cleanMeta),"ENSEMBL")]
print(dim(cleanData)) # 21940   365
print(dim(cleanMeta)) # 364   7
print(head(cleanMeta)) 
## save cleanData and cleanMeta
save(cleanData, cleanMeta, file = "Rdata/01.hsc.cleanData.cleanMeta.rdata" )


### ----------------------
### 03. Create Seurat obj
### ----------------------
load("Rdata/01.hsc.cleanData.cleanMeta.rdata")
dim(cleanData)
hsc <- CreateSeuratObject(counts = cleanData[,1:ncol(cleanData)-1], 
                                  meta.data = cleanMeta, 
                                  assay = "RNA", 
                                  min.cells = 3, 
                                  min.features = 200)
hsc <- PercentageFeatureSet(hsc, 
                           pattern = "^mt-", 
                           col.name = "pct_counts_Mito")
dim(hsc) # 11993 362
save(hsc,file = 'Rdata/temp01.hsc.seuratObj.rdata')
