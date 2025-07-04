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
setwd("/home/mayao/LABdata/MY/pdc")

### set Seurat
options(Seurat.object.assay.version = "v5")

## seurat save and load file
seurat_file <- 'Rdata/temp02.pdc.seurat.after_qc_norm.rdata'


### -----------------------
### 1. Filtering 
### -----------------------
##-----------------------##
load("Rdata/temp01.pdc.seuratObj.rdata")
## reset seurat obj name as 'obj' 
obj <- pdc
##-----------------------##
meta <- obj@meta.data %>% as.data.frame()

## filtering cutoff
low_dg <- 1500 
low_umi <- 2e+04
high_dg <- 1e+04
high_umi <- 3e+06
mito <- 3

## plot

##-----------------------##
group_var <- 'Name'
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
# HSC1  S12   ST 
#   26  116  132

obj <- subset(obj, subset=nCount_RNA > low_umi & 
                          nCount_RNA < high_umi &
                          nFeature_RNA > low_dg &
                          nFeature_RNA < high_dg &
                          pct_counts_Mito < mito)
dim(obj) # 11993   794
## After filtering, some genes might have rowSums = 0, because the only cells expressing these
## Check rowSums.
kept_genes <- list()
for (layer_n in Layers(obj[['RNA']])){
        counts <- obj$RNA[layer_n]
        genes <- counts[Matrix::rowSums(counts) > 0,] %>% rownames()
        kept_genes <- c(kept_genes,list(genes))
}

final_genes <- Reduce(intersect,kept_genes)
print(length(final_genes)) # 9681

for (layer_n in Layers(obj[['RNA']])){
        counts <- obj$RNA[layer_n]
        mx <- counts[final_genes,] 
        obj$RNA[layer_n] <- mx
}
obj <- obj[final_genes,]
print(dim(obj)) # 9681 794
table(obj$Name)
# Exp1 Exp2 Exp3 Exp4 Exp5 
#  144   90  379   27  154
# Exp1 Exp2 Exp3 Exp5 
#  144   90  406  154




# ### -----------------------
# ### 2. SCTransform 
# ### -----------------------
# obj <- SCTransform(obj, vars.to.regress='pct_counts_Mito')
# print(dim(obj))


# ### -----------------------
# ### 3. Cell cycle scoring
# ### -----------------------
# s.genes <- cc.genes$s.genes %>% str_to_title()
# g2m.genes <- cc.genes$g2m.genes %>% str_to_title()

# obj <- CellCycleScoring(obj, g2m.features=g2m.genes, s.features=s.genes)
# obj$cc_difference <- obj$S.Score - obj$G2M.Score
# obj$Phase <- factor(obj$Phase, levels=c("G1","S","G2M"))
# table(obj$Phase, obj$Group)
# #      ST S12 Fresh_HSC1
# #   G1   65  94         35
# #   S    65 100          5
# #   G2M  97  89          2

# #         Fresh_HSC1 HSC1 S12  ST
# #   G1          35  150  95  63
# #   S            5    2 132  84
# #   G2M          2    2 104 120


# ### --------------------------------------
# ### SCTransform after cellcycle scoring
# ### --------------------------------------
# obj <- SCTransform(obj, vars.to.regress=c('pct_counts_Mito', "cc_difference"))
# print(dim(obj)) # 10964 704


### --------------------------------------
### 4. Add isPair
### --------------------------------------
df <- obj@meta.data %>% as.data.frame()
df$isPair <- sapply(as.character(df$Pairs), function(x){
    if (x != 'Fresh'){
        if ( nrow(df[df$Pairs %in% x, ]) == 2 ){
            return('Paired')
        }else {
            return('Unpaired')
        }
    }else {
        return(x)
    }
})
df$isPair <- factor(df$isPair, levels = c('Paired', 'Unpaired', 'Fresh'))
obj$isPair <- df$isPair
obj@meta.data[obj$Group == "HSC1", ]$Group <- "Fresh_HSC1" 
obj$Group <- factor(obj$Group, levels=c("ST","S12", "Fresh_HSC1"))
table(obj$isPair, obj$Group)
#            ST S12 Fresh_HSC1
#   Paired   154 182          0
#   Unpaired  73 101          0
#   Fresh      0   0        194

#            ST S12 Fresh_HSC1
#   Paired   186 212          0
#   Unpaired  81 119          0
#   Fresh      0   0        196



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

obj[["RNA"]] <- split(obj[["RNA"]], f=obj$Name)
obj <- IntegrateLayers(
                    object = obj,
                    method = RPCAIntegration,
                    assay = "RNA",
                    orig.reduction = 'pca',
                    new.reduction = 'rpca',
                    k.weight = 20,
                    normalization.method="LogNormalize",
                    verbose=F,
                    dims = 1:8
                    )
obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
DimPlot(obj, reduction="rpca", group.by="Group")

## check normalize
df <- data.frame(col_sum=colSums(obj[["RNA"]]$data), Group=obj$Group)
my_compare <- list(c("Fresh_HSC1", "ST"), c("Fresh_HSC1", "S12"), c("ST", "S12"))
p2 <- ggplot(df, aes(x=Group, y=col_sum, color=Group))+
    geom_boxplot()+
    # scale_color_manual(values=colGroup)+
    stat_compare_means(comparisons=my_compare,label = "p.signif")+
    geom_quasirandom(aes(color=Group),size = 1, alpha =0.6)+
    stat_summary(fun.y=mean,aes(color=Group),geom='point',shape=5,size= 2)+
    theme_bw()+
    ylab("Aggragated expression")+
    ggtitle("")+
    force_panelsizes(rows=unit(7,"cm"), cols=unit(5,"cm"))
p2






##-----------------------##
pdc <- obj
save(pdc, file=seurat_file)
