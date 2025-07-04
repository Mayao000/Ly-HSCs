## conda env: monocle
rm(list=ls())
set.seed(1234)
options(warn = -1)

cat('*************************************\n')
dates <- as.character(Sys.Date())
cat('            ', dates, '             \n')
cat('*************************************\n')



### ------------------------------------------------
###                   Function Zone
### ------------------------------------------------
source("/home/mayao/script/Function.r")


### --------------------
### 01. Set up env
### --------------------
suppressMessages(library(dplyr))
library(ggplot2) %>% suppressMessages()
library(Seurat) %>% suppressMessages()
library(cowplot) %>% suppressMessages()
library(sctransform) %>% suppressMessages()
library(stringr) %>% suppressMessages()
library(Matrix) %>% suppressMessages()
library(reshape2) %>% suppressMessages()
library(scales) %>% suppressMessages()
library(ggrepel) %>% suppressMessages()
library(ggh4x)  %>% suppressMessages()

### set path
setwd("/home/mayao/LABdata/MY/hsc1_hsc2/mpp")

### set Seurat
options(Seurat.object.assay.version = "v5")

## seurat save and load file
seurat_file <- 'Rdata/temp02.mpp.seurat.after_qc_norm.rdata'


### -----------------------
### 1. Filtering 
### -----------------------
##-----------------------##
load("Rdata/temp01.mpp.seuratObj.rdata")
## reset seurat obj name as 'obj' 
obj <- mpp
##-----------------------##
meta <- obj@meta.data %>% as.data.frame()

## filtering cutoff
low_dg <- 1000 
low_umi <- 1e+04
high_dg <- 5000 
high_umi <- 5e+05
mito <- 3

## plot

##-----------------------##
group_var <- 'Group'
is.umi <- TRUE # Refers to raw umi counts

##-----------------------##
## 1. Distribution of counts per cell
p1<- meta %>% 
        ggplot(aes(x=nCount_RNA, fill=!!ensym(group_var), color=!!ensym(group_var))) + 
        geom_density(alpha = 0.2) + 
        theme_classic() +
        scale_x_log10() +
        geom_vline(xintercept=c(high_umi,low_umi))
p1

## 2. Distribution of detected genes per cell
p2<- meta %>% 
        ggplot(aes(x=nFeature_RNA, fill=!!ensym(group_var), color=!!ensym(group_var))) + 
        geom_density(alpha = 0.2) + 
        theme_classic() +
        scale_x_log10() +
        geom_vline(xintercept=c(high_dg,low_dg))
p2
## 3. Relationship between detected genes and counts
p3 <- FeatureScatter(obj, 
                    feature1 = "nCount_RNA", 
                    feature2 = "nFeature_RNA",
                    group.by=group_var)+
                    geom_vline(xintercept=c(high_umi,low_umi))+
                    geom_hline(yintercept=c(high_dg,low_dg))
p3
# 4. Relationship between mito.percent and counts
p4 <- FeatureScatter(obj, 
                    feature1 = "nCount_RNA",
                    feature2 = "pct_counts_Mito",
                    group.by=group_var)+
                    geom_vline(xintercept=c(high_umi,low_umi))+
                    geom_hline(yintercept=mito)
p4

##-----------------------##
cell_quality_pl <- c(list(p1),list(p2),list(p3),list(p4))
save_plotlist(cell_quality_pl, ncol=1, outfile=paste0("plot/QC.cutoff_check.cell_quality.pdf"), width=6, height=15)


##-----------------------##
## Filtering
removed <- subset(obj, subset=nCount_RNA < low_umi | nCount_RNA == low_umi |
                          nCount_RNA > high_umi | nCount_RNA == high_umi |
                          nFeature_RNA < low_dg | nFeature_RNA == low_dg |
                          nFeature_RNA > high_dg | nFeature_RNA == high_dg |
                          pct_counts_Mito > mito | pct_counts_Mito == mito)
message("#### Removed Cells' Group ####")
print(table(removed$Group))
# MPP5 MPP6 MPP7 
#    8    3    6

obj <- subset(obj, subset=nCount_RNA > low_umi & 
                          nCount_RNA < high_umi &
                          nFeature_RNA > low_dg &
                          nFeature_RNA < high_dg &
                          pct_counts_Mito < mito)
dim(obj) # 11686   239
## After filtering, some genes might have rowSums = 0, because the only cells expressing these
## Check rowSums.
table(rowSums(obj$RNA) == 0 ) # FALSE 11686
is.na(obj$CD135_median) %>% table()
# FALSE  TRUE 
#   250     6
obj$is_facs_null <- FALSE
obj@meta.data[is.na(obj$CD135_median), "is_facs_null"] <- TRUE 


### --------------------------------------
### LogNorm
### --------------------------------------
s.genes <- cc.genes$s.genes %>% str_to_title()
g2m.genes <- cc.genes$g2m.genes %>% str_to_title()

obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
obj <- NormalizeData(obj, assay="RNA")
obj <- FindVariableFeatures(obj, assay="RNA")

## cc
obj <- CellCycleScoring(obj, g2m.features=g2m.genes, s.features=s.genes)
obj$cc_difference <- obj$S.Score - obj$G2M.Score
obj$Phase <- factor(obj$Phase, levels=c("G1","S","G2M"))
table(obj$Phase, obj$Group)
#        MPP5 MPP6 MPP7
#   G1    41   48   42
#   S     24   11   28
#   G2M   13   22   10

## scale
obj <- ScaleData(obj, assay="RNA", vars.to.regress=c('pct_counts_Mito',"cc_difference"))

## PCA
obj <- RunPCA(obj, assay="RNA", npcs=30)
DimPlot(obj, reduction="pca", group.by="is_facs_null")
DimPlot(obj, reduction="pca", group.by="Group")



## check normalize
df <- data.frame(col_sum=colSums(obj[["RNA"]]$data), Group=obj$Group)
my_compare <- list(c("MPP5","MPP6"), c("MPP5","MPP7"), c("MPP6","MPP7"))
p2 <- ggplot(df, aes(x=Group, y=col_sum, color=Group))+
    geom_boxplot()+
    # scale_color_manual(values=colGroup)+
#     stat_compare_means(comparisons=my_compare,label = "p.signif")+
#     geom_quasirandom(aes(color=Group),size = 1, alpha =0.6)+
    stat_summary(fun.y=mean,aes(color=Group),geom='point',shape=5,size= 2)+
    theme_bw()+
    ylab("Aggragated expression")+
    ggtitle("")+
    force_panelsizes(rows=unit(7,"cm"), cols=unit(5,"cm"))
p2




##-----------------------##
mpp <- obj
# save(mpp, file=seurat_file)
save(mpp, file="Rdata/temp02.mpp.seurat.after_qc_norm.lognorm.rdata")
