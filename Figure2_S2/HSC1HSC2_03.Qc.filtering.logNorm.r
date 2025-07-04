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

### set path
setwd("/home/mayao/LABdata/MY/hsc1_hsc2/new")

### set Seurat
options(Seurat.object.assay.version = "v5")

## seurat save and load file
seurat_file <- 'Rdata/temp02.hsc.seurat.after_qc_norm.rdata'


### -----------------------
### 1. Filtering 
### -----------------------
##-----------------------##
load("Rdata/temp01.hsc.seuratObj.rdata")
## reset seurat obj name as 'obj' 
obj <- hsc
##-----------------------##
meta <- obj@meta.data %>% as.data.frame()

## filtering cutoff
low_dg <- 1000 
low_umi <- 1e+04
high_dg <- 4000 
high_umi <- 2e+05
mito <- 5

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
# HSC1 HSC2 
#   21    5

obj <- subset(obj, subset=nCount_RNA > low_umi & 
                          nCount_RNA < high_umi &
                          nFeature_RNA > low_dg &
                          nFeature_RNA < high_dg &
                          pct_counts_Mito < mito)
dim(obj) # 11993   336
## After filtering, some genes might have rowSums = 0, because the only cells expressing these
## Check rowSums.
table(rowSums(obj$RNA) == 0 ) # FALSE 11993



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

## scale
obj <- ScaleData(obj, assay="RNA", vars.to.regress=c('pct_counts_Mito',"cc_difference"))
obj <- RunPCA(obj, assay="RNA", npcs=30)
DimPlot(obj, reduction="pca", group.by="Group")

## check normalize
df <- data.frame(col_sum=colSums(obj[["RNA"]]$data), Group=obj$Group)
my_compare <- list(c("Fresh_HSC1", "ST"), c("Fresh_HSC1", "S12"), c("ST", "S12"))
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
hsc <- obj
# save(hsc, file=seurat_file)
save(hsc, file="Rdata/temp02.hsc.seurat.after_qc_norm.lognorm.rdata")
