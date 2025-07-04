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





# --------------------------------------------------------------------------------------------
## Loading hsc data
load("/home/mayao/LABdata/MY/hsc1_hsc2/new/Rdata/02.hsc.seuratobj.after_qc_normed.rdata")
# hsc <- hsc %>% subset(Group == "HSC1")
hsc$Name <- "hsc"
colnames(hsc) <- paste0("HSC.", colnames(hsc))
all(colnames(hsc) == rownames(hsc@meta.data))

load("Rdata/02.mpp.seuratobj.rdata")
mpp$Name <- "mpp"
colnames(mpp) <- paste0("MPP.", colnames(mpp))
all(colnames(mpp) == rownames(mpp@meta.data))

## Create combined Seurat object
counts_ls <- c(list(hsc[['RNA']]$counts), list(mpp[['RNA']]$counts))
names(counts_ls) <- c("hsc", "mpp")
meta_ls <- c(list(as.data.frame(hsc@meta.data[,c("Group",  "Name")])),
             list(as.data.frame(mpp@meta.data[,c("Group", "Name")])) )
meta <- Reduce(rbind, meta_ls)
combine <- CreateSeuratObject(counts=counts_ls,meta.data=meta)

## correct names
combine$Group <- gsub("MPP6", "N1", combine$Group)
combine$Group <- gsub("MPP5", "N2", combine$Group)
combine$Group <- gsub("MPP7", "Ly-I", combine$Group)
# combine$Group <- factor(combine$Group, levels = c("HSC1", "N1", "N2", "Ly-I"))
combine$Group <- factor(combine$Group, levels = c("HSC1", "HSC2", "N1", "N2", "Ly-I"))
table(combine$Group)
# HSC1   N1   N2 Ly-I  
#  159   81   78   80 
# HSC1 HSC2   N1   N2 Ly-I 
#  159  177   81   78   80 


## change colGroup
colGroup <- c("black", "red", 
              "#FF80C0","#188386", "#808080"
              )
colGroup <- setNames(colGroup, nm = levels(combine$Group))
colGroup

## calculate mito
combine <- PercentageFeatureSet(combine, 
			pattern = "^mt-", 
			col.name = "pct_counts_Mito")

## logNorm
DefaultAssay(combine) <- "RNA"
combine[["RNA"]] <- JoinLayers(combine[["RNA"]])
combine <- NormalizeData(combine, assay="RNA")

## cc
s.genes <- cc.genes$s.genes %>% str_to_title()
g2m.genes <- cc.genes$g2m.genes %>% str_to_title()
combine <- CellCycleScoring(combine, g2m.features=g2m.genes, s.features=s.genes)
combine$cc_difference <- combine$S.Score - combine$G2M.Score
combine$Phase <- factor(combine$Phase, levels=c("G1","S","G2M"))

## Set parameters after testing
nfeatures <- 1500
ndim <- 6
optRes <- 0.4

## Split the layers
combine[["RNA"]] <- split(combine[["RNA"]], f=combine$Name)

## Variable features 
combine <- FindVariableFeatures(combine, assay="RNA", nfeatures=nfeatures)

## Scale
combine <- ScaleData(combine, assay="RNA", vars.to.regress=c('pct_counts_Mito',"cc_difference"))

## PCA
combine <- RunPCA(combine, assay="RNA", npcs=30)
ElbowPlot(combine)
DimPlot(combine, group.by = "Group", reduction = "pca", pt.size = 1.5) +
      scale_color_manual(values = colGroup)


## Integragate
combine <- IntegrateLayers(
                        object = combine,
                        method = HarmonyIntegration,
                        assay = "RNA",
                        orig.reduction = 'pca',
                        new.reduction = "harmony",
                        # theta = ,  # Controls the diversity preservation. Higher values prioritize more diversity within batches. Typical values range from 1 to 10.,
                        # lamda = NULL, # Adjusts the importance of batch effects. A default value is often used. 
                        # sigma = 0.1, # Controls the bandwidth of Gaussian kernels in the algorithm. Often, the default is used.
                        # # normalization.method = normalization.method,
                        verbose=F,
                        dims = 1:ndim,
                        seed.use = 1234
                        )
DimPlot(combine, group.by = "Group", reduction = "harmony", pt.size = 1.5) +
      scale_color_manual(values = colGroup)

combine <- IntegrateLayers(
                        object = combine,
                        method = CCAIntegration,
                        assay = "RNA",
                        orig.reduction = 'pca',
                        new.reduction = "cca",
                        k.weight = 20, 
                        verbose=F,
                        dims = 1:5
                        )

combine <- IntegrateLayers(
                        object = combine,
                        method = RPCAIntegration,
                        assay = "RNA",
                        orig.reduction = 'pca',
                        new.reduction = "rpca",
                        k.weight = 20, 
                        verbose=F,
                        dims = 1:5
                        )

## Join the layers
combine[["RNA"]] <- JoinLayers(combine[["RNA"]])

obj <- combine
obj <- RunUMAP(obj, 
            reduction = "harmony", 
            reduction.name = "umap.harmony", 
            dims=1:6, 
            n.neighbors = 60,
            min.dist = 0.2,
            spread = 0.5,
            metric = "euclidean",
            seed.use = 1234,
            n.components = 3,
            verbose=FALSE)

library(plotly)
embed <- Embeddings(obj, reduction = "umap.harmony") %>% as.data.frame()
embed$Group <- obj$Group

p <- plot_ly(x=as.numeric(embed$umapharmony_1), y=embed$umapharmony_2, z=embed$umapharmony_3, type="scatter3d", mode="markers", color=embed$Group, colors = colGroup)
htmlwidgets::saveWidget(p, file = paste0("plot/96.RCombineHSC1.hsc1hsc2-n1n2n3.3D-umapharmony.html"))

embed <- embed %>% subset(Group != "HSC2")
p <- plot_ly(x=as.numeric(embed$umapharmony_1), y=embed$umapharmony_2, z=embed$umapharmony_3, type="scatter3d", mode="markers", color=embed$Group, colors = colGroup)
htmlwidgets::saveWidget(p, file = paste0("plot/96.RCombineHSC1.hsc1hsc2-n1n2n3.removeHSC2.3D-umapharmony.html"))

combine <- obj
save(combine, file = "Rdata/RCombineHSC1.02.hsc-n1n2n3.combine-obj.rdata")


## Plot
# load("Rdata/RCombineHSC1.02.hsc-n1n2n3.combine-obj.rdata")
reduction_name <- "umap.harmony"
reduction_1 <- paste0(gsub(".", "", reduction_name, fixed=T),"_1")
reduction_2 <- paste0(gsub(".", "", reduction_name, fixed=T),"_2")

umap <- FetchData(combine, vars = c(reduction_1, reduction_2, "seurat_clusters", "Group"))
p1 <- ggplot(umap, aes(x=!!ensym(reduction_1), y=!!ensym(reduction_2), fill=Group))+
        geom_point(size = 1.5, alpha=1, shape = 21, color = "white", stroke = 0.2) +
      #   geom_point(size = 1, alpha=1) +
        scale_fill_manual(values = colGroup, name="") +
        plotTheme +
        labs(x="UMAP 1", y="UMAP 2") +
        guides(color = guide_legend(override.aes = list(size=1.5))) +
        ggtitle("") +
        theme(axis.text = element_text(color = "black", size = 12), 
              axis.title = element_text(color = "black", size = 12))+
        force_panelsizes(rows=unit(4,"cm"), cols=unit(5,"cm"))
p1


umap <- FetchData(combine %>% subset(Group != "HSC2"), vars = c(reduction_1, reduction_2, "seurat_clusters", "Group"))
p2<- ggplot(umap, aes(x=!!ensym(reduction_1), y=!!ensym(reduction_2), fill=Group))+
        geom_point(size = 1.5, alpha=1, shape = 21, color = "white", stroke = 0.2) +
      #   geom_point(size = 1, alpha=1) +
        scale_fill_manual(values = colGroup, name="") +
        plotTheme +
        labs(x="UMAP 1", y="UMAP 2") +
        guides(color = guide_legend(override.aes = list(size=1.5))) +
        ggtitle("") +
        theme(axis.text = element_text(color = "black", size = 12), 
              axis.title = element_text(color = "black", size = 12))+
        force_panelsizes(rows=unit(4,"cm"), cols=unit(5,"cm"))
p2


pdf("plot/89.hsc1hsc2-n1n2n3.umap.pdf", width = 10, height = 10)
print(plot_grid(plotlist=c(list(p1), list(p2)), ncol = 2))
dev.off()

