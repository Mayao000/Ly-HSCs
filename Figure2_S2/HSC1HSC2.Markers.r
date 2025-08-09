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


### ----------------------------------------------
### Find Markers
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





