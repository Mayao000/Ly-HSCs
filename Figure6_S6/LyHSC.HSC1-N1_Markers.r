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

### set plot options
### Define color palettes and plot themes
colLibrary <- colorRampPalette(brewer.pal(n = 5, name = "Spectral"))(5)
colGEX <- c("grey85", brewer.pal(7, "Reds"))      # color for gene expression
colCcy <- c("black", "blue", "darkorange")        # color for cellcycle phase
plotTheme <- theme_classic(base_size = 14)
colGroup <- c("black", 
              "#FF80C0","#188386", "#808080"
              )
show_col(colGroup)

## set other options after being decided
colCls <- c("#4189C9", "#D791BF","#FFCF3F")
show_col(colCls)
colState <- c("3" = "#CE95C2", "2" = "#4F67B0", "7" = "#07A6ED", "1" = "#74D109",
              "5" = "#D1272A", "6" = "#FBAA19")
show_col(colState)

## set python options and load modules
python_path <- '/home/mayao/.conda/envs/seurat/bin/'
use_python(python_path)
py_run_string("import numpy as np")
py_run_string("import matplotlib.pyplot as pl")

### set path
setwd("/home/mayao/LABdata/MY/hsc1_hsc2/mpp/")

### set Seurat
options(Seurat.object.assay.version = "v5")



## ------------------------------------------------------------------------------------
# load("Rdata/RCombineHSC1.05.hsc-n1n2n3.combine-obj.monocle-info-added.rdata")
load("Rdata/RCombineHSC1.02.hsc-n1n2n3.combine-obj.rdata")
obj <- combine

## ------------------------------------------------------------------------------------
## Markers between HSC1 and HSC2-N1
group.by <- "Group"
Idents(obj) <- group.by 

res_ls <- my_markers(obj, find.all=FALSE, ident.1 = "N1", ident.2 = "HSC1", group.by = group.by, whether_plot = FALSE)
markers <- res_ls[["markers"]]
# markers_clean <- res_ls[["markers_clean"]]
# go <- res_ls[["enrich_res"]][["go"]]
# kegg <- res_ls[["enrich_res"]][["kegg"]]

# xlsx_file <- "20.1.hsc-n1n2n3.N1-HSC1.markers.xlsx"
# xlsx_file <- "20.2.hsc-n1n2n3.N1-HSC1.markers-full.xlsx"
xlsx_file <- "table/20.3.hsc1hsc2-n1n2n3.N1-HSC1.markers-full.xlsx"
if ( file.exists(paste0("table/",xlsx_file)) ) {
  file.remove( paste0("table/",xlsx_file) )
}
wb <- my_xlsx_create(xlsx_file)

wb <- my_xlsx_write(wb, data = markers, sheet = "markers", rowNames = TRUE)
# wb <- my_xlsx_write(wb, data = markers_clean, sheet = "markers", rowNames = TRUE)
# wb <- my_xlsx_write(wb, data = go, sheet = "GO", rowNames = TRUE)
# wb <- my_xlsx_write(wb, data = kegg, sheet = "KEGG", rowNames = TRUE)
my_xlsx_save(wb, xlsx_file)

View(markers)
nrow(markers)
markers$Direction <- ifelse(markers$avg_log2FC > 0, "HSC2-N1", "HSC1")
table(markers$Direction)




## -----------------------------------------------------------------------------------
## box plot
library(ggpubr)
library(ggbeeswarm)
library(RColorBrewer)
sub_obj <- obj %>% subset(Group %in% c("HSC1", "N1"))

features_n1 <- c("Myl10",  "Plac8",  "Ltb",  "Tmsb10",  "Ier2",  "Ramp1",  "Cd47",  "Arpc1b",  "Ifitm1",  "Actb",  "Tmsb4x",  "Ifitm2",  "Cfl1",  "Ptma")
features_hsc1 <- c("Serinc3",  "Npm1",  "Sord",  "Mmrn1",  "Acadl",  "Upp1",  "Itgb3",  "Tgm2",  "Apoe",  "Trim47",  "Itga2b",  "Clec1a",  "Aldh1a1",  "Sult1a1")

features <- c(features_n1, rev(features_hsc1))
boxdata <- FetchData(sub_obj, vars=c(features,"Group"))
boxdata_mut <- tidyr::pivot_longer(boxdata, cols = -Group, names_to = "Feature", values_to = "Expression")
boxdata_mut$Feature <- factor(boxdata_mut$Feature, levels = features)
table(boxdata_mut$Feature)

my_comparisons <- list(c("HSC1", "N1"))
p <- ggplot(boxdata_mut, aes(x=Group, y=Expression, fill=Group))+
        geom_boxplot(color = "black", outlier.shape = NA, alpha = 1) + 
        facet_wrap( vars(Feature), ncol = 7 )+
        scale_fill_manual(values = c("darkgrey", "#FF80C0"))+
        stat_compare_means(comparisons=my_comparisons, 
                           method="wilcox.test",
                           label = "p.signif",
                           label.y = max(boxdata_mut$Expression) - 1)+
        labs(y="Expression",x="")+ 
        scale_color_manual(values = c("darkgrey", "#FF80C0"))+
        theme_bw() + 
        theme(panel.border = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.line = element_line(colour = "black"),
                legend.position="none", 
                axis.text = element_text(color = "black", size = 20),
                strip.text = element_text(color = "black", size = 20, face = "italic"))+
        force_panelsizes(rows=unit(3,"cm"), cols=unit(3,"cm"))
ggsave(paste0("plot/87.hsc-n1n2n3.markers.N1-HSC1.boxplot.pdf"), width = 15, height = 15)
