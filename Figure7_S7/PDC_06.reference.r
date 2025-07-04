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
            #   "#E02C18",
              "navajowhite3", "darkorange1","forestgreen","royalblue1",
              "#188386", "#FF80C0", "#808080"
              )
colGroup <- setNames(colGroup, nm = c("HSC1", paste0("MPP", seq(7))))
show_col(colGroup)

## set other options after being decided
colCls <- pal_npg("nrc")(10)[seq(7)]
show_col(pal_npg("nrc")(10))


## set python options and load modules
python_path <- '/home/mayao/.conda/envs/seurat/bin/'
use_python(python_path)
py_run_string("import numpy as np")
py_run_string("import matplotlib.pyplot as pl")

### set path
setwd("/home/mayao/LABdata/MY/hsc1_hsc2/mpp/")

### set Seurat
options(Seurat.object.assay.version = "v5")



# ## ----------------------------------------------------------------------------------------
# ## Loading hsc data
# load("/home/mayao/LABdata/MY/hsc1_hsc2/new/Rdata/02.hsc.seurat.clustered.markered.rdata")
# count_hsc <- hsc@assays[["RNA"]]$counts
# meta_hsc <- hsc@meta.data[,c("Group","Marker")] %>% as.data.frame()
# meta_hsc$Name <- "HSC"
# meta_hsc$CD135_median <- "Null"
# meta_hsc$CD34_median <- "Null"
# colnames(count_hsc) <- paste0("HSC.", colnames(count_hsc))
# rownames(meta_hsc) <- paste0("HSC.", rownames(meta_hsc))
# all(colnames(count_hsc) == rownames(meta_hsc))

# ##---------------------------------------------------------------------
# ## loading DF data
# load("/home/mayao/LABdata/DF/01.UMI.TPM.Meta.Homeostasis.Rdata") # count TPM meta.bg
# dim(TPM) # 28579  1270
# dim(count) # 28579  1270
# dim(meta.bg) # 1270   11
# table(meta.bg$Sorting.cell)
# meta.bg$Marker <- meta.bg$Marker %>% as.vector()
# meta.bg[meta.bg$Sorting.cell == "ST_HSC", ]$Marker %>% table()
# # CD34-Flt3-Lin-Kit+Sca1+ CD34+Flt3+Lin-Kit+Sca1+ 
# #                      22                      19 
# meta.bg[meta.bg$Sorting.cell == "MPP1", ]$Marker %>% table()
# # CD34-Flt3-Lin-Kit+Sca1+CD150+CD48- 
# #                                72
# meta.bg[meta.bg$Sorting.cell == "MPP2", ]$Marker %>% table()
# # CD34-Flt3-Lin-Kit+Sca1+CD150+CD48+ 
# #                                 72
# meta.bg[meta.bg$Sorting.cell == "MPP3", ]$Marker %>% table()
# # CD34-Flt3-Lin-Kit+Sca1+CD150-CD48+ 
# #                                 76
# meta.bg[meta.bg$Sorting.cell == "MPP4", ]$Marker %>% table()
# # CD34+Flt3+Lin-Kit+Sca1+CD150-CD48+ 
# #                                 78
# meta.bg[meta.bg$Sorting.cell == "LT_HSC", ]$Marker %>% table()
# # CD34-Flt3-Lin-Kit+Sca1+ 
# #                      22

# ## Extract MPP1,2,3,4,ST_HSC,LT_HSC,LMPP
# celltypes <- levels(meta.bg$Sorting.cell)
# kepttypes <- celltypes[-which(celltypes %in% c("ST_HSC","LT_HSC","LMPP"))]
# cells_df <- meta.bg %>% 
#     # subset(subset=Sorting.cell %in% c("MPP1", "MPP2", "MPP3", "MPP4", "LT_HSC", "ST_HSC", "LMPP")) %>% 
#     subset(subset=Sorting.cell %in% kepttypes) %>% 
#     rownames()
# length(cells_df ) # 1176
# count_df <- count[,cells_df ] %>% as.matrix()
# dim(count_df) # 28579   1176
# meta_df <- meta.bg[cells_df, ]
# dim(meta_df) # 1176  11
# meta_df <- meta_df[, c("Sorting.cell", "Marker")]
# colnames(meta_df)[1] <- "Group"
# meta_df$Name <- "DF"
# meta_df$CD135_median <- "Null"
# meta_df$CD34_median <- "Null"


# ## load MPP6/5/7
# load("Rdata/02.mpp.seuratobj.rdata")
# table(mpp$Group)
# # MPP6 MPP5 MPP7 
# #   81   78   80
# count_mpp <- mpp@assays[["RNA"]]$counts 
# meta_mpp <- mpp@meta.data[,c("Group","Marker", "CD135_median", "CD34_median")] %>% as.data.frame()
# meta_mpp$Name <- "MPP"
# colnames(count_mpp) <- paste0("MPP.", colnames(count_mpp))
# rownames(meta_mpp) <- paste0("MPP.", rownames(meta_mpp))
# all(colnames(count_mpp) == rownames(meta_mpp))


# ## keep only common genes
# all(rownames(count_df) == rownames(count_hsc) )
# all(rownames(count_hsc) == rownames(count_mpp) )
# temp <- intersect(rownames(count_df) , rownames(count_hsc) ) ## 11928
# kept_genes <- intersect(temp, rownames(count_mpp)) ## 11162
# count_df <- count_df[kept_genes, ]
# count_hsc <- count_hsc[kept_genes,]
# count_mpp <- count_mpp[kept_genes,]

# ## combine
# count_ls <- c(list(count_df), list(count_hsc), list(count_mpp))
# names(count_ls) <- c("DF", "HSC", "MPP")
# meta <- rbind(meta_df, meta_hsc, meta_mpp)
# meta$Group <- factor(meta$Group, levels=c("HSC1", "Fraction I", "Fraction II", "Fraction III",
#                                           "ESLAM", "ESLAMSK",
#                                           "HSC2","MPP6", "MPP5", "MPP7",
#                                           "MPP1","MPP2", "MPP3","MPP4",
#                                           "HPC2", "HPC3",
#                                           "CMP", "CLP","GMP", "MEP",
#                                           "MK", "EryA", "EryB", 
#                                           "B cell", "CD4T", "CD8T", "NK cell",
#                                           "Granulocyte", "Monocyte", "Macrophage"))  
# table(meta$Group)
# #         HSC1   Fraction I  Fraction II Fraction III        ESLAM      ESLAMSK 
# #          159          102           49           66           58           57 
# #         HSC2         MPP6         MPP5         MPP7         MPP1         MPP2 
# #          177           81           78           80           72           72 
# #         MPP3         MPP4         HPC2         HPC3          CMP          CLP 
# #           76           78           44           45           22           30 
# #          GMP          MEP           MK         EryA         EryB       B cell 
# #           24           19           19           37           33           40 
# #         CD4T         CD8T      NK cell  Granulocyte     Monocyte   Macrophage 
# #           46           42           39           41           41           24 



# ## --------------------------------------------------------------------------------------------------
# ## 2. Seurat obj
# blood_new <- CreateSeuratObject(counts = count_ls, meta.data = meta, min.cells = 0, min.features = 0)
# table(blood_new$Group)
# save(blood_new, file="Rdata/temp01.blood_new.before_any_process.rdata")
# print(dim(blood_new)) # 11162  1751


## --------------------------------------------------------------------------------------------------
## 3. Process only HSC & MPP
load("Rdata/temp13.blood_new.before_any_process.notfiltergenes.rdata")
obj <- blood_new %>% subset(subset = Group %in% c("HSC1",  paste0("MPP", seq(7))))
obj$Group <- factor(obj$Group, levels = c("HSC1", paste0("MPP", seq(7))))
table(obj$Group)
#  HSC1 MPP1 MPP2 MPP3 MPP4 MPP5 MPP6 MPP7 
#  159   72   72   76   78   78   81   80


## pct_counts_Mito
obj <- PercentageFeatureSet(obj, 
			pattern = "^mt-", 
			col.name = "pct_counts_Mito")

## log norm
s.genes <- cc.genes$s.genes %>% str_to_title()
g2m.genes <- cc.genes$g2m.genes %>% str_to_title()
obj <- NormalizeData(obj, assay="RNA")

# ## Variable features
obj <- FindVariableFeatures(obj, assay="RNA", nfeature = 3000)

## cc
obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
obj <- CellCycleScoring(obj, g2m.features=g2m.genes, s.features=s.genes)
obj$cc_difference <- obj$S.Score - obj$G2M.Score
obj$Phase <- factor(obj$Phase, levels=c("G1","S","G2M"))
table(obj$Phase, obj$Group)
#       HSC1 HSC2 MPP1 MPP2 MPP3 MPP4 MPP5 MPP6 MPP7
#   G1   155  176   30    3    0    5   59   78   58
#   S      3    1   29   36   33   27   18    3   21
#   G2M    1    0   13   33   43   46    1    0    1


## scale
obj <- ScaleData(obj, assay="RNA", vars.to.regress=c('pct_counts_Mito',"cc_difference"))

# ## PCA
obj <- RunPCA(obj, assay="RNA", npcs=30)
DimPlot(obj, reduction="pca", group.by="Group", pt.size = 1.5)+
    scale_color_manual(values = colGroup)
ElbowPlot(obj) # 5-8

## UMAP
obj <- RunUMAP(obj, 
            reduction = "pca", 
            reduction.name = "umap", 
            dims=1:6, 
            verbose=FALSE)
DimPlot(obj, reduction="umap", group.by="Group", pt.size = 1.5)+
    scale_color_manual(values = colGroup)


## --------------------------------------------------------------------------------------------------
## 3. Process only HSC & MPP
## Decide parameters
assay <- "RNA"
normalization.method <- "LogNormalize"
nfeature <- 1500
ndim <- 5

## load data
load("Rdata/temp13.blood_new.before_any_process.notfiltergenes.rdata")
obj <- blood_new %>% subset(subset = Group %in% c("HSC1",  paste0("MPP", seq(7))))
obj$Group <- factor(obj$Group, levels = c("HSC1", "MPP6", "MPP5", "MPP7", paste0("MPP", seq(4))))
table(obj$Group)

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

## Split
obj[[assay]] <- split(obj[[assay]], f=obj$Name)

## Variable features
obj <- FindVariableFeatures(obj, assay="RNA", nfeature = nfeature)

# ## scale
obj <- ScaleData(obj, assay="RNA", vars.to.regress=c('pct_counts_Mito',"cc_difference"))

hspc <- obj 
save(hspc, file = "Rdata/03.combine.hspc.seuratobj.rdata")




