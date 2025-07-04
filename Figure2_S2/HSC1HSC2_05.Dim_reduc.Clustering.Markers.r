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
colGroup <- c("black", "#E02C18")
show_col(colGroup)
colTissue <- brewer.pal(4,'RdGy')[c(1,4)]
colCondition <- brewer.pal(4,'PuOr')[c(3,4)]

## set other options after being decided
nPC <- 5
nClust <- 11
colCls <- c("#CAB2D6", "#6A3D9A", "grey") 

## set python options and load modules
python_path <- '/home/mayao/.conda/envs/seurat/bin/'
use_python(python_path)
# py_run_string("import numpy as np")
# py_run_string("import matplotlib.pyplot as pl")



### set path
setwd("/home/mayao/LABdata/MY/hsc1_hsc2/new")

### set Seurat
options(Seurat.object.assay.version = "v5")


### -----------------------
### 1. Dimension reduction
### -----------------------
## Reduc: PCA, t-SNE, UMAP

##-----------------------##
load('Rdata/temp02.hsc.seurat.after_qc_norm.lognorm.rdata')
obj <- hsc
dim(hsc)
table(hsc$Group)
# HSC1 HSC2 
#  159  177

##-----------------------##
## PCA
obj <- RunPCA(obj, npcs = 30, verbose = F)
## Plot PCA
reduc <- 'pca'
p1 <- dimplot(obj, reduction = reduc, group.by = 'Group')
p2 <- dimplot(obj, reduction = reduc, group.by = 'Condition')
p3 <- dimplot(obj, reduction = reduc, group.by = 'Phase')
p4 <- dimplot(obj, reduction = reduc, group.by = 'Tissue')
p5 <- ElbowPlot(obj) # nPC set to 5
pca_pl <- c(list(p1), list(p2), list(p3),list(p4),list(p5))
## DimHeatmap
library(viridis)
dimheatmap <- DimHeatmap(hsc, dims = 1, balanced = T, nfeatures = 30, fast = F) + 
                scale_fill_viridis()+
                force_panelsizes(rows = unit(3, "cm"), cols = unit(3, "cm"))+
                theme(axis.text = element_text(size = 16, face = 'italic'), 
                    legend.text = element_text(size = 16))
dimheatmap                    

##-----------------------##
## UMAP
obj <- RunUMAP(obj, dims = 1:nPC, verbose = F)
## Plot UMAP
reduc <- 'umap'
p1 <- dimplot(obj, reduction = reduc, group.by = 'Group')
p2 <- dimplot(obj, reduction = reduc, group.by = 'Condition')
p3 <- dimplot(obj, reduction = reduc, group.by = 'Phase')
p4 <- dimplot(obj, reduction = reduc, group.by = 'Tissue')
umap_pl <- c(list(p1), list(p2), list(p3),list(p4))

##-----------------------##
## tSNE
obj <- RunTSNE(obj, dims = 1:nPC, verbose = F)
## Plot tSNE
reduc <- 'tsne'
p1 <- dimplot(obj, reduction = reduc, group.by = 'Group')
p2 <- dimplot(obj, reduction = reduc, group.by = 'Condition')
p3 <- dimplot(obj, reduction = reduc, group.by = 'Phase')
p4 <- dimplot(obj, reduction = reduc, group.by = 'Tissue')
tsne_pl <- c(list(p1), list(p2), list(p3),list(p4))

##-----------------------##
out <- c(pca_pl,umap_pl,tsne_pl)

pdf(paste0("plot/01.hsc.Dim_reduc.no_cls.pdf"), width=15, height=42)
print(plot_grid(pca_pl[[1]], pca_pl[[2]], pca_pl[[3]],
                pca_pl[[4]], pca_pl[[5]], NULL,
                umap_pl[[1]], umap_pl[[2]], umap_pl[[3]],
                umap_pl[[4]], NULL, NULL,
                tsne_pl[[1]], tsne_pl[[2]], tsne_pl[[3]],
                tsne_pl[[4]], NULL, NULL,
                ncol=3, byrow=TRUE))
dev.off()

##-----------------------## 
## save for ai
ggsave(filename="plot/z01.hsc.grp.pca.pdf", plot=pca_pl[[1]])
ggsave(filename="plot/z01.hsc.grp.umap.pdf", plot=umap_pl[[1]])
ggsave(filename="plot/z01.hsc.grp.tsne.pdf", plot=tsne_pl[[1]])
ggsave(filename="plot/z01.hsc.grp.pca.dimheatmap.pdf", plot=dimheatmap)


## save
hsc <- obj
table(hsc$Group)
save(hsc, file = "Rdata/02.hsc.seuratobj.after_qc_normed.rdata")



### -----------------------
### 3. Find Markers
### -----------------------
## Group
obj$Group <- factor(obj$Group, levels = c("HSC1","HSC2"))
Idents(obj) <- "Group"
markers <- FindAllMarkers(obj, only.pos=TRUE, verbose=FALSE) %>% 
    subset(subset=p_val_adj < 0.05) %>%
    plotly::arrange(cluster, desc(avg_log2FC))
View(markers)
obj@misc$markers_grp <- markers
table(markers$cluster)

enrich_res <- my_enrich(markers, directed=FALSE, universe=rownames(obj),pvalueCutoff=0.05, qvalueCutoff=1, go_cat="BP")
go <- enrich_res[["go"]]
kegg <- enrich_res[["kegg"]]
View(go)
View(kegg)

## save results
xlsx_file <- '02.hsc.markers.grp.250225-new.xlsx'
wb <- my_xlsx_create(xlsx_file)
wb <- my_xlsx_write(wb, sheetName="Group_markers", data=markers, rowNames=TRUE)
wb <- my_xlsx_write(wb, sheetName="GO", data=go, rowNames=TRUE)
my_xlsx_save(wb, xlsx_file)





