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
library(ggsci) %>% suppressMessages()
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
library(ggalluvial) %>% suppressMessages()
library(clusterProfiler) %>% suppressMessages()
library(cols4all)  %>% suppressMessages()

### set plot options
### Define color palettes and plot themes
colLibrary <- colorRampPalette(brewer.pal(n = 5, name = "Spectral"))(5)
colGEX <- c("grey85", brewer.pal(7, "Reds"))      # color for gene expression
colCcy <- c("black", "blue", "darkorange")        # color for cellcycle phase
plotTheme <- theme_classic(base_size = 14)
colGroup <- c("black", 
            #   "#E02C18",
              "navajowhite3", "darkorange1","forestgreen","royalblue1",
              "#188386", "#FF80C0", "#808080"
              )
show_col(colGroup)
colGroup <- setNames(colGroup, nm = c("HSC1", paste0("MPP", seq(7))))

colRef <- c("black", 
              "navajowhite3", "darkorange1","forestgreen","royalblue1",
              "#188386", "#FF80C0", "#808080"
              )
colRef <- setNames(colRef, nm = c("HSC1", paste0("MPP", seq(7))))

## set other options after being decided
colCls <- pal_npg("nrc")(10)[c(1,2,5,4,3)]

## set python options and load modules
python_path <- '/home/mayao/.conda/envs/seurat/bin/'
use_python(python_path)
py_run_string("import numpy as np")
py_run_string("import matplotlib.pyplot as pl")

### set path
setwd("/home/mayao/LABdata/MY/hsc1_hsc2/mpp/")

### set Seurat
options(Seurat.object.assay.version = "v5")

## predict label order
predict_label_order <- c("HSC1", paste0("MPP", seq(7)))



### ---------------------------------------------------------------------------------------
### Data and genesets
### ---------------------------------------------------------------------------------------
load("/home/mayao/genome_file/GeneSets/Konturek.gs_ls.rdata")
load("/home/mayao/genome_file/GeneSets/2017_NAR.2020_CSC.2021_Blood.gene_set_list.rdata")
names(gene_set_list)
names(gs_ls)

load("Rdata/02.mpp.seuratobj.rdata")
obj <- mpp
Idents(obj) <- "Group"

### ---------------------------------------------------------------------------------------
### Scoring
### ---------------------------------------------------------------------------------------
seuratobj <- obj %>% subset(subset = Group %in% c("MPP5", "MPP6","MPP7"))
gene_set_list <- gs_ls
geneset_name <- names(gs_ls)
group.by <- "Group"
assay <- "RNA"
meta_list <- list()
for (id in geneset_name){
    message("#### ",id,"####")
    genes_to_plot <- gene_set_list[[id]]
    genes_to_plot <- genes_to_plot[genes_to_plot %in% rownames(seuratobj)]
    exp_mat <- as.matrix(seuratobj[[assay]]$data[genes_to_plot,])
    print(dim(exp_mat))
    meta <- seuratobj@meta.data %>% plotly::select(!!ensym(group.by))
    print(dim(meta))
    meta <- bind_cols(meta, as.data.frame(t(exp_mat)))
    if (ncol(meta) > 1){
    meta <- tidyr::pivot_longer(meta, -group.by, names_to="Gene", values_to="Expression")
    meta_summary <- meta %>%
                    group_by(!!ensym(group.by), Gene) %>%
                    plotly::summarise(Avg = mean(Expression))
    meta_aggr <- meta_summary %>% 
                    group_by(!!ensym(group.by)) %>%
                    plotly::summarise(Avg_aggr = sum(Avg)/length(genes_to_plot), 
                        Feature=id)
    meta_list[[id]] <- meta_aggr
    }
}
df <- Reduce(rbind, meta_list)
df$Feature <- factor(df$Feature, levels=unique(df$Feature))


### ---------------------------------------------------------------------------------------
### Radar
library(ggradar)
data_wide <- tidyr::pivot_wider(df, names_from = Feature, values_from = Avg_aggr, values_fill = 0)
data_wide
data_scaled <- data_wide %>%
  mutate(across(-Group, ~ (. - min(.)) / (max(.) - min(.)))) 

ggradar(data_scaled, 
        base.size = 8,
        group.line.width = 1, 
        group.point.size = 0.2, 
        legend.position = "bottom", 
        fill = T, 
        fill.alpha = 0.2) +
    scale_color_manual(values = colGroup) +
    scale_fill_manual(values = colGroup) +
    labs(title = "")
ggsave("plot/43.mpp.Radar.Konturek")
### ---------------------------------------------------------------------------------------




