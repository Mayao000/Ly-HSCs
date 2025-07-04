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
library(plyr) %>% suppressMessages()
library(cowplot) %>% suppressMessages()
library(sctransform) %>% suppressMessages()
library(stringr) %>% suppressMessages()
library(clustree) %>% suppressMessages()
library(RColorBrewer) %>% suppressMessages()
library(scales) %>% suppressMessages()
library(patchwork)  %>% suppressMessages()
library(ggh4x)  %>% suppressMessages()
library(Matrix) %>% suppressMessages()
library(ggrepel) %>% suppressMessages()
library(data.table) %>% suppressMessages()
library(reticulate) %>% suppressMessages()

### set plot options
### Define color palettes and plot themes
colLibrary <- colorRampPalette(brewer.pal(n = 5, name = "Spectral"))(5)
colGEX <- c("grey85", brewer.pal(7, "Reds"))      # color for gene expression
colCcy <- c("black", "blue", "darkorange")        # color for cellcycle phase
plotTheme <- theme_classic(base_size = 14)
# colGroup <- colorRampPalette(brewer.pal(n = 4, name = "Spectral"))(2)
colGroup <- c("#FF80C0", "#188386", "#808080")
show_col(colGroup)
colTissue <- brewer.pal(4,'RdGy')[c(1,4)]
colCondition <- brewer.pal(4,'PuOr')[c(3,4)]

## set other options after being decided
colCls <- colorRampPalette(brewer.pal(n = 3, name = "Set1"))(3)

## set python options and load modules
python_path <- '/home/mayao/.conda/envs/seurat/bin/'
use_python(python_path)
# py_run_string("import numpy as np")
# py_run_string("import matplotlib.pyplot as pl")



### set path
setwd("/home/mayao/LABdata/MY/hsc1_hsc2/mpp")

### set Seurat
options(Seurat.object.assay.version = "v5")






### ---------------------------------------------------------------------------------
### 1. Dimension reduction and clustering
load('Rdata/temp02.mpp.seurat.after_qc_norm.lognorm.rdata')
obj <- mpp
obj$Group <- factor(obj$Group, levels = c("MPP6", "MPP5", "MPP7"))

## some settings 
## other_var, color and genesets
other_var <- "Group"
colGroup_new <- setNames(colGroup, nm=levels(obj$Group))
colList <- c(list(colGroup_new))
names(colList) <- "Group"
load('/home/mayao/genome_file/GeneSets/2017_NAR.2020_CSC.2021_Blood.gene_set_list.rdata')
names(gene_set_list)


### Decided parameters
assay <- "RNA"
vars_to_regress <- c('pct_counts_Mito', "cc_difference")
nfeature <- 2500
ndim <- 6
method <- "pca"
res <- 0.5

## Scale and PCA
obj <- FindVariableFeatures(obj, assay=assay, nfeature=nfeature, verbose=FALSE)
obj <- ScaleData(obj, assay=assay, vars.to.regress=vars_to_regress, verbose=FALSE)
obj <- RunPCA(obj, assay=assay, npcs=30, verbose=FALSE)

### Clustering
obj <- FindNeighbors(obj, reduction = method, dims = 1:ndim, verbose = FALSE)
obj <- FindClusters(obj, resolution = res, algorithm = 4, verbose = FALSE)

## UMAP
reduction_name <- paste0("umap")
obj <- RunUMAP(obj, 
            reduction = method, 
            reduction.name = reduction_name, 
            dims=1:ndim, 
            verbose=FALSE)

## Clusters
obj$seurat_clusters <- NULL
obj$seurat_clusters <- obj[[paste0('RNA_snn_res.', as.character(res))]]
obj$seurat_clusters <- factor(obj$seurat_clusters)
Idents(obj) <- obj$seurat_clusters
table(obj$seurat_clusters)
# 1  2  3 
# 94 92 53


## save
mpp <- obj
save(mpp, file = "Rdata/02.mpp.seuratobj.rdata")


## Plot
## Plot final umap plot
umap <- FetchData(obj, vars = c("umap_1", "umap_2", "seurat_clusters", "Group"))
colVar_ls <- c(list(colCls),list(colGroup))
names(colVar_ls) <- c("seurat_clusters", "Group")
p_ls <- list()
for (var in c("seurat_clusters", "Group")){
    colVar <- colVar_ls[[var]]
    p_ls[[var]] <- ggplot(umap, aes(umap_1, umap_2, color=!!ensym(var)))+
                    geom_point(size = 2) +
                    scale_color_manual(values = colVar, name="") +
                    plotTheme +
                    labs(x="UMAP 1", y="UMAP 2") +
                    guides(color = guide_legend(override.aes = list(size=3))) +
                    force_panelsizes(rows=unit(8,"cm"), cols=unit(8,"cm"))
}
pdf(paste0("plot/01.mpp.nfeature-2500.ndim-6.res-05.umap.pdf"), width = 20, height = 20)
print(plot_grid(plotlist=p_ls, ncol = 2))
dev.off()
