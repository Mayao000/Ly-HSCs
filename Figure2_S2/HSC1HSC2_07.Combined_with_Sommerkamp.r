## conda env: monocle
rm(list=ls())
set.seed(1234)
options(warn = -1, future.globals.maxSize= 2*1024^3)

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

## set python options and load modules
python_path <- '/home/mayao/.conda/envs/seurat/bin/'
use_python(python_path)
py_run_string("import numpy as np")
py_run_string("import matplotlib.pyplot as pl")
## set other opts
skPC <- 6
skClust <- 6
# colCls <- colorRampPalette(brewer.pal(n = 10, name = "Paired"))(skClust)
# colCls[4] <- "brown4"
# show_col(colCls)
colCls <- c("#D23F27", "#D0B38B", "#AC6FAE", 
            "#93A0D1", "#A4D071", "#F9A442")
colCcy <- c("#DC6D73", "#E5A0A0", "#8C7CB9", "#DFC5E0")

# colorRampPalette(brewer.pal(n = 10, name = "Paired"))(10) %>% show_col()
# colorRampPalette(brewer.pal(n = 10, name = "Paired"))(10)

## set plot opts
plotTheme <- theme_classic(base_size = 14)
colGroup <- c("black", "#E02C18")
colGroup_new <- c(colGroup, colorRampPalette(brewer.pal(n = 10, name = "Paired"))(10))
show_col(colGroup_new)

## setwd
setwd("/home/mayao/LABdata/MY/hsc1_hsc2/new")




## ************************************************************************************
## load data
load("Rdata/02.hsc.seuratobj.after_qc_normed.rdata")
load("Rdata/Pub01.Sommerkamp.sk.cleaned.seurat.rdata")
hsc$Name <- "hsc"
sk$Name <- "sk"


## Create combined Seurat object
counts_ls <- c(list(hsc[['RNA']]$counts), list(sk[['RNA']]$counts))
names(counts_ls) <- c("hsc", "sk")
meta_ls <- c(list(as.data.frame(hsc@meta.data[,c("Group",  "Name")])), 
             list(as.data.frame(sk@meta.data[,c("Group",  "Name")])))
meta <- Reduce(rbind, meta_ls)
combine <- CreateSeuratObject(counts=counts_ls,meta.data=meta)


## calculate mito
combine <- PercentageFeatureSet(combine, 
			pattern = "^mt-", 
			col.name = "pct_counts_Mito")

## logNorm
DefaultAssay(combine) <- "RNA"
s.genes <- cc.genes$s.genes %>% str_to_title()
g2m.genes <- cc.genes$g2m.genes %>% str_to_title()

combine[["RNA"]] <- JoinLayers(combine[["RNA"]])
combine <- NormalizeData(combine, assay="RNA")

## cc
combine <- CellCycleScoring(combine, g2m.features=g2m.genes, s.features=s.genes)
combine$cc_difference <- combine$S.Score - combine$G2M.Score
combine$Phase <- factor(combine$Phase, levels=c("G1","S","G2M"))

## color and genesets
load('/home/mayao/genome_file/GeneSets/2017_NAR.2020_CSC.2021_Blood.gene_set_list.rdata')


## set parameters after testing
nfeatures <- 2000
ndim <- 6
integrate_method <- "cca"
wt <- 25
# colGroup_new <- setNames(colGroup_new, nm=levels(combine$Group_new))

# ## Variable features 
combine <- FindVariableFeatures(combine, assay="RNA", nfeatures=nfeatures)

# ## scale
combine <- ScaleData(combine, assay="RNA", vars.to.regress=c('pct_counts_Mito',"cc_difference"))

## PCA
combine <- RunPCA(combine, assay="RNA", npcs=30)

# Integrate
combine[["RNA"]] <- split(combine[["RNA"]], f=combine$Name)
combine <- IntegrateLayers(
                        object = combine,
                        method = CCAIntegration,
                        assay = "RNA",
                        orig.reduction = 'pca',
                        new.reduction = integrate_method,
                        k.weight = wt,
                        normalization.method="LogNormalize",
                        verbose=F,
                        dims = 1:ndim
                        )
combine[["RNA"]] <- JoinLayers(combine[["RNA"]])

### Clustering
optRes <- 0.4
combine <- FindNeighbors(combine, reduction=integrate_method, dims = 1:ndim, verbose = FALSE)
combine <- FindClusters(combine, resolution = optRes, algorithm = 4, verbose = FALSE)

## UMAP
reduction_name <- paste0("umap.", integrate_method)
combine <- RunUMAP(combine, 
            reduction = integrate_method, 
            reduction.name = reduction_name, 
            dims=1:ndim, 
            verbose=FALSE)

## Clusters
combine$seurat_clusters <- NULL
combine$seurat_clusters <- combine[[paste0('RNA_snn_res.', as.character(optRes))]]
combine$seurat_clusters <- factor(combine$seurat_clusters)
Idents(combine) <- combine$seurat_clusters
table(combine$seurat_clusters)

DimPlot(combine, reduc = "umap.cca")

## re-name clusters
meta <- combine@meta.data %>% as.data.frame()
combine$Cluster_new <- sapply(as.character(colnames(combine)), function(x){
        if (meta[x, "seurat_clusters"] == "5"){
                Cluster_new <- "C1-HSC"
        }else if (meta[x, "seurat_clusters"] == "4"){
                Cluster_new <- "C2-MPP1"
        }else if (meta[x, "seurat_clusters"] == "2"){
                Cluster_new <- "C3-MPP5"
        }else if (meta[x, "seurat_clusters"] == "3"){
                Cluster_new <- "C6-MPP2/3"
        }else if (meta[x, "seurat_clusters"] == "1"){
                Cluster_new <- "C5-MPP3"
        }else if (meta[x, "seurat_clusters"] == "6"){
                Cluster_new <- "C4-MPP4/5"
        }
        return(Cluster_new)
})
combine$Cluster_new <- factor(combine$Cluster_new)
table(combine$Cluster_new)


## save
save(combine, file="Rdata/03.hsc-sk.combine.rdata")


## Plot
load("Rdata/03.hsc-sk.combine.rdata")
reduction_name <- "umap.cca"
reduction_1 <- paste0(gsub(".", "", reduction_name, fixed=T),"_1")
reduction_2 <- paste0(gsub(".", "", reduction_name, fixed=T),"_2")
group.by.new <- "Cluster_new"
highlight_by <- "Group"
highlight_grp <- c("HSC1","HSC2")
colHighlight <- setNames(colGroup,nm = c("HSC1","HSC2"))
umap <- FetchData(combine, vars = c(reduction_1, reduction_2, group.by.new, highlight_by))
subdata <- umap[ umap[,highlight_by] == highlight_grp ,]
p1 <- ggplot(umap, aes(x=!!ensym(reduction_1), y=!!ensym(reduction_2), color=!!ensym(group.by.new) ))+
        geom_point(stroke = 0.4, size = 0.1, alpha=0.2) +
        geom_point(data=subdata, aes(x=!!ensym(reduction_1), y=!!ensym(reduction_2), fill=!!ensym(highlight_by)),shape=21, color="white", stroke = 0.5, size = 3) +
        scale_color_manual(values = colCls, name="") +
        scale_fill_manual(values = as.character(colHighlight[highlight_grp]), name="", label=highlight_grp) +
        plotTheme +
        labs(x="UMAP 1", y="UMAP 2") +
        guides(color = guide_legend(override.aes = list(size=3))) +
        ggtitle("") +
        force_panelsizes(rows=unit(8,"cm"), cols=unit(9,"cm"))
p1


p2 <- ggplot(umap, aes(x=!!ensym(reduction_1), y=!!ensym(reduction_2), color=!!ensym(group.by.new)))+
        geom_point(size = 0.1, alpha=0.9) +
        scale_color_manual(values = colCls, name="") +
        plotTheme +
        labs(x="UMAP 1", y="UMAP 2") +
        guides(color = guide_legend(override.aes = list(size=3))) +
        ggtitle("") +
        force_panelsizes(rows=unit(8,"cm"), cols=unit(9,"cm"))
p2
pdf("plot/27.hsc-combine.umap-featureplot.pdf", width = 30, height = 10)
print(plot_grid(plotlist=c(list(p1), list(p2)), ncol = 3))
dev.off()


load('/home/mayao/genome_file/GeneSets/2017_NAR.2020_CSC.2021_Blood.gene_set_list.rdata')
group.by.new <- "Cluster_new"
colors <- c("#222A8A", "#D3CAB1", "#FFA83E", "#B70335")
geneset_name <- c("HSC", paste0("MPP", seq(1,5)))
module_res <- geneset_dotplot(seuratobj=combine, assay="RNA", gene_set_list=gene_set_list, group.by = group.by.new,
                                geneset_name = geneset_name, row = "feature", max_by = "Group", 
                                highlight_method="score")
p_module <- module_res[['p']] + 
                scale_fill_gradientn(colors = colors)+
                scale_x_discrete(labels=c("C1", "C2", "C3", "C4", "C5", "C6"))+
                theme(legend.position="bottom", legend.direction="vertical")+
                force_panelsizes(rows=unit(6.5,"cm"), cols=unit(7,"cm"))
p_module
ggsave("plot/27.hsc-combine.p_module.cluster.featureplot.pdf")


### Markers Plots
load("Rdata/03.hsc-sk.combine.rdata")
### FACS markers
features <- c("Kit","Ly6a","Procr", "Slamf1","Cd48","Cd34","Flt3")
group.by <- "Cluster_new"
p <- DotPlot(combine, features = features, group.by=group.by) +
        theme_bw() +
        # scale_color_gradientn(colors=c("black","black","red")) +
        scale_color_distiller(palette = "RdGy")+
        scale_size(range=c(2,8))+
        labs(x="",y="")+
        coord_flip() +
        scale_y_discrete(labels=c("C1", "C2", "C3", "C4", "C5", "C6"))+
        theme(axis.text.x=element_text(angle=45, hjust=1))+
        # theme(legend.position = "bottom", legend.direction = "horizontal")+
        force_panelsizes(rows=unit(9,"cm"), cols=unit(6,"cm"))
p
ggsave(p, width = 10, height = 10, filename = "plot/29.hsc-combine.facs-markers.dotplot.pdf")

### Lineage markers --- dotplot
features <- c("Vwf","Gata1", "Klf1", "Epor", "Elane", "Cebpe", "Ctsg", "Mpo", "Gfi1", "Dntt", "Cd52")
group.by <- "Cluster_new"
p <- DotPlot(combine, features = features, group.by=group.by) +
        theme_bw() +
        # scale_color_gradientn(colors=c("grey","#E41A1C")) +
        scale_color_distiller(palette = "RdGy")+
        scale_size(range=c(2,8))+
        labs(x="",y="")+
        coord_flip() +
        scale_y_discrete(labels=c("C1", "C2", "C3", "C4", "C5", "C6"))+
        theme(axis.text.x=element_text(angle=45, hjust=1))+
        force_panelsizes(rows=unit(8,"cm"), cols=unit(7,"cm"))
p
ggsave(p, width = 10, height = 10, filename = "plot/30.hsc-combine.lineage-markers.dotplot.pdf")