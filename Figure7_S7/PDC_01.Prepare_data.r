# '''
# This is to analyse pdc.
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
setwd("/home/mayao/LABdata/MY/pdc")


# ### --------------------
# ### 1. Prepare data
# ### --------------------
## 1.1 process meta
# meta <- read.table("meta/PDC_meta_new.txt", header = T)
# meta$Library <- sapply(meta$Library, function(x){
#     gsub('lib','Library',x)
# })
# meta$Cell <- sapply(meta$Cell, function(x){
#     gsub('lib','Library',x)
# })
# rownames(meta) <- meta$Cell
# meta$Group <- factor(meta$Group, levels = c('ST', 'S12','Fresh_HSC1'))
# meta$Batch <- factor(meta$Batch)
# meta$Pairs <- factor(meta$Pairs)
# meta$Condition <- factor(meta$Condition, levels = c('Cultured', 'Fresh'))
# write.table(meta, file='meta/PDC_meta_new_revised.txt', sep = '\t', row.names = F, quote = F)

## 
meta <- read.table('meta/PDC_meta_new_revised.txt', header=TRUE)
rownames(meta) <- meta$Cell
meta_1 <- meta %>% subset(subset = Library %in% paste0("Library",seq(1,5)))
meta_2 <- meta %>% subset(subset = Library %in% paste0("Library",seq(6,11)))
meta_3 <- meta %>% subset(subset = Library %in% paste0("Library",seq(12,13)))

## 1.2 remove non-coding and duplicated genes
count_1 <- read.table("count/PDC_1_5.rawCount.txt", header = T)
count_1 <- count_1[, c("ENSEMBL", "geneTypes", "genes", rownames(meta_1))]
cleanData_1 <- keepCodingGenes(count_1) %>% removeDuplicatedGenes() 

count_2 <- read.table("count/PDC_6_11.rawCount.txt", header = T)
count_2 <- count_2[, c("ENSEMBL", "geneTypes", "genes", rownames(meta_2))]
cleanData_2 <- keepCodingGenes(count_2) %>% removeDuplicatedGenes()

count_3 <- read.table("count/PDC_12_13.rawCount.txt", header = T)
count_3 <- count_3[, c("ENSEMBL", "geneTypes", "genes", rownames(meta_3))]
cleanData_3 <- keepCodingGenes(count_3) %>% removeDuplicatedGenes()

cleanMeta_1 <- meta_1[colnames(cleanData_1[1:ncol(cleanData_1)-1]) ,] 
cleanMeta_1 <- cleanMeta_1[!is.na(cleanMeta_1$Group),]
cleanData_1 <- cleanData_1[, c(rownames(cleanMeta_1),"ENSEMBL")]

cleanMeta_2 <- meta_2[colnames(cleanData_2[1:ncol(cleanData_2)-1]) ,] 
cleanMeta_2 <- cleanMeta_2[!is.na(cleanMeta_2$Group),]
cleanData_2 <- cleanData_2[, c(rownames(cleanMeta_2),"ENSEMBL")]

cleanMeta_3 <- meta_3[colnames(cleanData_3[1:ncol(cleanData_3)-1]) ,] 
cleanMeta_3 <- cleanMeta_3[!is.na(cleanMeta_3$Group),]
cleanData_3 <- cleanData_3[, c(rownames(cleanMeta_3),"ENSEMBL")]

print(dim(cleanData_1)) # pdc data 1-5: 21940  434
print(dim(cleanMeta_1)) # pdc meta 1-5: 433   10

print(dim(cleanData_2)) # pdc data 6-11: 21940   523
print(dim(cleanMeta_2)) # pdc meta 6-11: 522  10

print(dim(cleanData_3)) # pdc data 12-13: 21940   151
print(dim(cleanMeta_3)) # 150 10


## 1.3 save cleanData and cleanMeta
save(cleanData_1, cleanMeta_1,
     cleanData_2, cleanMeta_2, 
     cleanData_3, cleanMeta_3,
     file = "Rdata/01.pdc.cleanData.cleanMeta.rdata" )


### ----------------------
### 03. Create Seurat obj
### ----------------------
## load pdc 
load("Rdata/01.pdc.cleanData.cleanMeta.rdata")
exp1 <- cleanData_1[,grepl("Library1|Library2", colnames(cleanData_1))] %>% as.matrix()
colnames(exp1) <- paste0("pdc.", colnames(exp1))
exp2 <- cleanData_1[,grepl("Library3|Library4|Library5", colnames(cleanData_1))] %>% as.matrix()
colnames(exp2) <- paste0("pdc.", colnames(exp2))
exp3 <- cleanData_2[,1:ncol(cleanData_2)-1] %>% as.matrix()
colnames(exp3) <- paste0("pdc.", colnames(exp3))
exp4 <- cleanData_3[,1:ncol(cleanData_3)-1] %>% as.matrix()
colnames(exp4) <- paste0("pdc.", colnames(exp4))

rownames(cleanMeta_1) <- paste0("pdc.", rownames(cleanMeta_1))
rownames(cleanMeta_2) <- paste0("pdc.", rownames(cleanMeta_2))
rownames(cleanMeta_3) <- paste0("pdc.", rownames(cleanMeta_3))

## load hsc1
load("/home/mayao/LABdata/MY/hsc1_hsc2/new/Rdata/01.hsc.cleanData.cleanMeta.rdata")
hsc1 <- cleanMeta[cleanMeta$Group == "HSC1", ] %>% rownames()
exp5 <- cleanData[,hsc1] %>% as.matrix()
colnames(exp5) <- paste0("hsc1.", colnames(exp5))
cleanMeta <- cleanMeta[hsc1,]
rownames(cleanMeta) <- paste0("hsc1.", rownames(cleanMeta))

colnames(cleanMeta)[grepl("Marker",colnames(cleanMeta))] <- "Sorting"
cleanMeta$Pairs <- "Fresh"
cleanMeta$Date <- ""
cleanMeta$Batch <- 7

count_ls <- c(list(exp1),list(exp2),list(cbind(exp3,exp4)), list(exp5))
names(count_ls) <- c("Exp1","Exp2","Exp3","Exp5")

meta_ls <- c(list(cleanMeta_1), list(cleanMeta_2), list(cleanMeta_3), list(cleanMeta))
meta <- Reduce(rbind,meta_ls)
meta$Library <- str_split(rownames(meta), "[_]", simplify=TRUE)[,1]

meta$Name <- sapply(as.character(meta$Library), function(x){
    if (grepl("Library6|Library7|Library8|Library9|Library10|Library11", x)){
        return("Exp3")
    }else if (grepl("pdc.Library3|pdc.Library4|pdc.Library5", x)){
        return("Exp2")
    }else if (x == "pdc.Library1" | x == "pdc.Library2"){
        return("Exp1")
    }else if (grepl("pdc.Library12|pdc.Library13", x)){
        return("Exp3")
    }else{
        return("Exp5")
    }
})
table(meta$Name)
# Exp1 Exp2 Exp3 Exp4 Exp5
#  165  268  522  150  182
## after merge exp3 and exp4
# Exp1 Exp2 Exp3 Exp5 
#  165  268  672  182


pdc <- CreateSeuratObject(counts = count_ls, 
                                  meta.data = meta, 
                                  assay = "RNA", 
                                  min.cells = 3, 
                                  min.features = 200)
pdc <- PercentageFeatureSet(pdc, 
                           pattern = "^mt-", 
                           col.name = "pct_counts_Mito")
dim(pdc) # 13712  1068
table(pdc$Name)
# Exp1 Exp2 Exp3 Exp4 Exp5 
#  160  140  505   83  180
## after merging exp3 and exp4
# Exp1 Exp2 Exp3 Exp5 
#  160  140  588  180
save(pdc,file = 'Rdata/temp01.pdc.seuratObj.rdata')


