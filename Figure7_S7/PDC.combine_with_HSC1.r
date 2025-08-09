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
colLibrary <- colorRampPalette(brewer.pal(n = 11, name = "Spectral"))(11)
colGEX <- c("grey85", brewer.pal(7, "Reds"))      # color for gene expression
colCcy <- c("black", "blue", "darkorange")        # color for cellcycle phase
plotTheme <- theme_classic(base_size = 14)
# colGroup <- colorRampPalette(brewer.pal(n = 4, name = "Spectral"))(2)
colGroup <- c("#2171B5","#ED85B0", "#A1D99B")
show_col(colGroup)

colTissue <- brewer.pal(4,'RdGy')[c(1,4)]
colCondition <- brewer.pal(4,'PuOr')[c(3,4)]
show_col(colCondition)

colBatch <- colorRampPalette(brewer.pal(n = 7, name = "Set1"))(7)
show_col(colBatch)

colisPair <- c("#74C476", "#BAE4B3", "#9E9AC8")
show_col(colisPair)

colMode <- c("steelblue1", "plum1", "lightgrey", "#74C476")
show_col(colMode)

colName <- c("red", "blue", "green", "purple")


## set other options after being decided
nPC <- 10
nClust <- 3
# colCls <- colorRampPalette(brewer.pal(n = nClust, name = "Spectral"))(nClust)
colCls <- c("#52A052", "#B79EFF", "#FA0087", "#3264B0")
show_col(colCls)

## set python options and load modules
python_path <- '/home/mayao/.conda/envs/seurat/bin/'
use_python(python_path)
# py_run_string("import numpy as np")
# py_run_string("import matplotlib.pyplot as pl")



### set path
setwd("/home/mayao/LABdata/MY/pdc")

### set Seurat
options(Seurat.object.assay.version = "v5")

## seurat save file
seurat_save <- 'Rdata/temp03.pdc.seurat.dim_reduc.clustered.rdata'


### -----------------------
### 1. Dimension reduction
### -----------------------
## Reduc: PCA, t-SNE, UMAP

##------------------------------------------------------------------------------
load('Rdata/temp02.pdc.seurat.after_qc_norm.rdata')
obj <- pdc

## pct_counts_Mito
obj <- PercentageFeatureSet(obj, 
			pattern = "^mt-", 
			col.name = "pct_counts_Mito")

## log norm
s.genes <- cc.genes$s.genes %>% str_to_title()
g2m.genes <- cc.genes$g2m.genes %>% str_to_title()
obj <- NormalizeData(obj, assay="RNA")

## cc
obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
obj <- CellCycleScoring(obj, g2m.features=g2m.genes, s.features=s.genes)
obj$cc_difference <- obj$S.Score - obj$G2M.Score
obj$Phase <- factor(obj$Phase, levels=c("G1","S","G2M"))
table(obj$Phase, obj$Group)



## Set parameters after testing
nfeatures <- 1000
ndim <- 4
integrate_method <- "cca"
wt <- 25
optRes <- 0.3

## Split the layers
obj[["RNA"]] <- split(obj[["RNA"]], f=obj$Name)

## Variable features 
obj <- FindVariableFeatures(obj, assay="RNA", nfeatures=nfeatures)

## Scale
obj <- ScaleData(obj, assay="RNA", vars.to.regress=c('pct_counts_Mito',"cc_difference"))

## PCA
obj <- RunPCA(obj, assay="RNA", npcs=30)

# Integrate
obj <- IntegrateLayers(
                        object = obj,
                        method = CCAIntegration,
                        assay = "RNA",
                        orig.reduction = 'pca',
                        new.reduction = integrate_method,
                        k.weight = wt,
                        normalization.method="LogNormalize",
                        verbose=F,
                        dims = 1:ndim
                        )

## Join the layers
obj[["RNA"]] <- JoinLayers(obj[["RNA"]])

### Clustering
obj <- FindNeighbors(obj, reduction=integrate_method, dims = 1:ndim, verbose = FALSE)
obj <- FindClusters(obj, resolution = optRes, algorithm = 4, verbose = FALSE)

## UMAP
reduction_name <- paste0("umap.", integrate_method)
obj <- RunUMAP(obj, 
            reduction = integrate_method, 
            reduction.name = reduction_name, 
            dims=1:ndim, 
            verbose=FALSE)

## Clusters
obj$seurat_clusters <- NULL
obj$seurat_clusters <- obj[[paste0('RNA_snn_res.', as.character(optRes))]]
obj$seurat_clusters <- factor(obj$seurat_clusters)
DimPlot(obj, reduction = reduction_name, group.by = "seurat_clusters")/
DimPlot(obj, reduction = reduction_name, group.by = "Group")

## re-name clusters
meta <- obj@meta.data %>% as.data.frame()
seurat_clusters <- sapply(as.character(colnames(obj)), function(x){
        if (meta[x, "seurat_clusters"] == "2"){
                return("C1") 
        }else if (meta[x, "seurat_clusters"] == "4"){
                return("C2") 
        }else if (meta[x, "seurat_clusters"] == "1"){
                return("C3") 
        }else if (meta[x, "seurat_clusters"] == "3"){
                return("C4") 
        }
})
obj$seurat_clusters <- seurat_clusters
obj$seurat_clusters <- factor(obj$seurat_clusters)
table(obj$seurat_clusters)

Idents(obj) <- obj$seurat_clusters
table(obj$seurat_clusters)

pdc <- obj
save(pdc, file = "Rdata/03.pdc.withHSC1.clustered.rdata")
##------------------------------------------------------------------------------




##------------------------------------------------------------------------------
## Plot
load("Rdata/03.pdc.withHSC1.clustered.rdata")
obj <- pdc
reduction_name <- "umap.cca"
reduction_1 <- paste0(gsub(".", "", reduction_name, fixed=T),"_1")
reduction_2 <- paste0(gsub(".", "", reduction_name, fixed=T),"_2")

umap <- FetchData(obj, vars = c(reduction_1, reduction_2, "seurat_clusters", "Group", "RF_"))
p1 <- ggplot(umap, aes(x=!!ensym(reduction_1), y=!!ensym(reduction_2), fill=Group))+
        geom_point(size = 2, alpha=1, shape = 21, color = "white", stroke = 0.5) +
        scale_fill_manual(values = colGroup, name="") +
        plotTheme +
        labs(x="UMAP 1", y="UMAP 2") +
        guides(color = guide_legend(override.aes = list(size=1.5))) +
        ggtitle("") +
        theme(axis.text = element_text(color = "black", size = 12), 
              axis.title = element_text(color = "black", size = 12))+
        force_panelsizes(rows=unit(5,"cm"), cols=unit(6,"cm"))
p1
p2 <- ggplot(umap, aes(x=!!ensym(reduction_1), y=!!ensym(reduction_2), fill=seurat_clusters))+
        geom_point(size = 2, alpha=1, shape = 21, color = "white", stroke = 0.5) +
        scale_fill_manual(values = colCls, name="") +
        plotTheme +
        labs(x="UMAP 1", y="UMAP 2") +
        guides(color = guide_legend(override.aes = list(size=1.5))) +
        ggtitle("") +
        theme(axis.text = element_text(color = "black", size = 12), 
              axis.title = element_text(color = "black", size = 12))+
        force_panelsizes(rows=unit(5,"cm"), cols=unit(6,"cm"))
p2
# ### plot components
library(gtools)
library(ggalluvial)
sub <- obj
orig_df <- sub@meta.data %>% as.data.frame()
table(orig_df$Group, orig_df$seurat_clusters)

df_count <- plyr::ddply(orig_df, .(seurat_clusters, Group), summarize, Counts = length(Group)) 
df_perc <- plyr::ddply(df_count, .(seurat_clusters), summarize,  Group = Group, Percent = Counts/sum(Counts)*100) 
df_perc$seurat_clusters <- factor(df_count$seurat_clusters, levels = levels(df_count$seurat_clusters))
table(df_perc$seurat_clusters)

p3 <- ggplot(df_perc, aes(x = seurat_clusters, y = Percent, fill = Group, stratum = Group, alluvium = Group)) + 
        geom_flow(width = 0.5, alpha = 0.3, knot.pos=0, color = 'white') + 
        geom_col(width = 0.5, color = 'white', alpha = 0.9) + 
        scale_y_continuous(expand = c(0, 0)) +  
        scale_fill_manual(values = colGroup) +
        xlab("") + 
        ylab("Explained variation (%)") + 
        theme_classic() +
        theme(legend.title = element_blank(),
              axis.text = element_text(color = "black", size = 12), 
              axis.title = element_text(color = "black", size = 12))+
        force_panelsizes(rows=unit(5,"cm"), cols=unit(5,"cm"))
p3

pdf("plot/71.pdc.umap-cca.1000-4-25-03.components.pdf", width = 30, height = 10)
# pdf("plot/temp.pdf", width = 30, height = 10)
print(plot_grid(plotlist=c(list(p1), list(p2), list(p3)), ncol = 3))
dev.off()
##------------------------------------------------------------------------------
