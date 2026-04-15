rm(list=ls())
set.seed(1234)
source("/home/mayao/script/Function.r")

## ***************************************************************************************************************
### Set up env
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

### set path
setwd("/home/mayao/LABdata/MY/hsc1_hsc2/new")

## set colors
colGroup <- c("black", "#E02C18")



## ***************************************************************************************************************
## HSC1 and HSC2
## ***************************************************************************************************************
############ Pre-rocess
meta <- read.table("meta/HSC1_HSC2_my.meta.txt", header = T)
rownames(meta) <- meta$Cell
origData <- read.table("matrix/HSC1_HSC2_my.rawCount.txt", header = T)
cleanData <- keepCodingGenes(origData) %>% removeDuplicatedGenes()  ## only coding genes and no duplication
                                                                    # column names: lib1_sc01 lib1_sc02 ... lib8_sc96 ENSEMBL
cleanMeta <- meta[colnames(cleanData[1:ncol(cleanData)-1]) ,] 
cleanMeta <- cleanMeta[!is.na(cleanMeta$Group),]
print(cleanMeta[1:30,])
cleanData <- cleanData[, c(rownames(cleanMeta),"ENSEMBL")]
print(dim(cleanData)) # 21940   365
print(dim(cleanMeta)) # 364   7
print(head(cleanMeta)) 
## save cleanData and cleanMeta
save(cleanData, cleanMeta, file = "Rdata/01.hsc.cleanData.cleanMeta.rdata" )


############ Create Seurat obj
load("Rdata/01.hsc.cleanData.cleanMeta.rdata")
dim(cleanData)
hsc <- CreateSeuratObject(counts = cleanData[,1:ncol(cleanData)-1], 
                                  meta.data = cleanMeta, 
                                  assay = "RNA", 
                                  min.cells = 3, 
                                  min.features = 200)
hsc <- PercentageFeatureSet(hsc, 
                           pattern = "^mt-", 
                           col.name = "pct_counts_Mito")
dim(hsc) # 11993 362
save(hsc,file = 'Rdata/temp01.hsc.seuratObj.rdata')


############ QC 
load("Rdata/temp01.hsc.seuratObj.rdata")
obj <- hsc
meta <- obj@meta.data %>% as.data.frame()

## filtering cutoff
low_dg <- 1000 
low_umi <- 1e+04
high_dg <- 4000 
high_umi <- 2e+05
mito <- 5

removed <- subset(obj, subset=nCount_RNA < low_umi | nCount_RNA == low_umi |
                          nCount_RNA > high_umi | nCount_RNA == high_umi |
                          nFeature_RNA < low_dg | nFeature_RNA == low_dg |
                          nFeature_RNA > high_dg | nFeature_RNA == high_dg |
                          pct_counts_Mito > mito | pct_counts_Mito == mito)
obj <- subset(obj, subset=nCount_RNA > low_umi & 
                          nCount_RNA < high_umi &
                          nFeature_RNA > low_dg &
                          nFeature_RNA < high_dg &
                          pct_counts_Mito < mito)


############ LogNorm
s.genes <- cc.genes$s.genes %>% str_to_title()
g2m.genes <- cc.genes$g2m.genes %>% str_to_title()

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

hsc <- obj
save(hsc, file="Rdata/temp02.hsc.seurat.after_qc_norm.lognorm.rdata")


############ Dimension reduction
load('Rdata/temp02.hsc.seurat.after_qc_norm.lognorm.rdata')
obj <- hsc
## UMAP
ElbowPlot(obj) # nPC set to 5
obj <- RunUMAP(obj, dims = 1:nPC, verbose = F)
hsc <- obj
save(hsc, file = "Rdata/02.hsc.seuratobj.after_qc_normed.rdata")

## Plot umap
## *********** Figure 2A ***********
umap_df <- FetchData(obj, vars=c("umap_1", "umap_2", "Group"))
ggplot(umap_df, aes(x=umap_1, y=umap_2, color=Group)) +
    geom_point(size=1.2)+
    scale_color_manual(values=colGroup,name="") +
    theme_classic() +
    labs(x="UMAP 1", y="UMAP 2") +
    theme(text = element_text(color="black")) +
    force_panelsizes(rows=unit(5,"cm"), cols=unit(5,"cm"))
ggsave("plot/z06.hsc.umap.pdf")



############ Markers
## *********** Figure 2B ***********
library(ComplexHeatmap)
library(circlize)
## Group
obj$Group <- factor(obj$Group, levels = c("HSC1","HSC2"))
Idents(obj) <- "Group"
markers <- FindAllMarkers(obj, only.pos=TRUE, verbose=FALSE) %>% 
    subset(subset=p_val_adj < 0.05) %>%
    plotly::arrange(cluster, desc(avg_log2FC))
markers_clean <- markers[!grepl("Rps|Rpl|mt",markers$gene),] # do not show ribosomal and mito genes
obj@misc$markers_clean <- markers_clean

od <- order(obj$Group)
mx <- obj[["RNA"]]$data[markers_clean$gene, od]  # 470 598
mx <- zScore(mx)

cols <- setNames(colGroup[1:2], nm=c("HSC1", "HSC2"))
cell_annotation <- HeatmapAnnotation(df = data.frame(row.names=colnames(mx),
                                              Group=obj@meta.data[,"Group"][od]),
                                    col = list(Group=cols),
                                    which = "col")

breaks <- c(-1, 0, 1)  # Custom break points
hp_col <- colorRamp2(breaks,c("#73A8D9", "white", "#F7CF2B"))

hsc1_gene <- c("Cavin2", "Slamf1", "Nupr1", "Selp", "Sult1a1", "Esam", "Mllt3","Gata2", "Pdzk1ip1")
hsc2_gene <- c("Tespa1", "Ccl3", "Flt3", "Plac8", "Cd52", "Cd53", "Cd34", "Il17ra", "Cd27")
hsc1_gene_index <- match(hsc1_gene, rownames(mx))
hsc2_gene_index <- match(hsc2_gene, rownames(mx))
gene_index <- c(hsc1_gene_index,hsc2_gene_index)
genes <- c(hsc1_gene, hsc2_gene)
mark_colors <- rep(colGroup, c(length(hsc1_gene), length(hsc2_gene)))
right_annotation <- rowAnnotation(foo = anno_mark(at = gene_index, labels = genes, 
                                                     labels_gp = gpar(fontface = "italic", col=mark_colors)),
                                          width = unit(6, "cm"))

p <- Heatmap(mx, 
    col=hp_col,
    name="Z-Score",
    top_annotation=cell_annotation,
    right_annotation = right_annotation,
    cluster_rows=TRUE, 
    cluster_columns=FALSE, 
    show_row_names=FALSE,
    show_column_names=FALSE,
    show_row_dend=FALSE,
    heatmap_legend_param = list(direction = "horizontal"),
    heatmap_width = unit(12, "cm"),
    heatmap_height = unit(10, "cm"))
p
pdf("plot/z09.hsc.markers.heatmap.pdf", width=15, height=15)
print(p)
dev.off()



############ GO and KEGG
enrich_res <- my_enrich(markers, 
                        directed = FALSE,
                        log2fc_cutoff = 0,
                        universe=rownames(obj),  
                        pvalueCutoff=0.05, 
                        qvalueCutoff=1,
                        GO = TRUE,
                        go_cat = c("BP"),
                        KEGG = TRUE)
go <- enrich_res[['go']] %>% as.data.frame()
kegg <- enrich_res[['kegg']] %>% as.data.frame()

## *********** Figure 2C ***********
selected_term_index <- c(1,4,7,19,20,21,24,26,29)
col_bar <- setNames(c("lightskyblue", "lightcoral"), nm = c("HSC1", "HSC2"))
p_list <- my_enrich_barplot(plot_df = go, 
                            directed = FALSE,
                            cluster_level = c("HSC1", "HSC2"), 
                            log2fc_cutoff = 0,
                            selected_term_index = selected_term_index,
                            colors = col_bar,
                            bar_width = 0.5,
                            x_var="qvalue",
                            y_var="Description")
p1 <- p_list[[1]] + 
        ylab("")+
        xlab("-log10(q-value)")+
        ggtitle("HSC1") +
        theme(text = element_text(color="black",family="sans"),
            axis.text.x = element_text(size=8))+
        force_panelsizes(row=unit(2.5,"cm"),col=unit(8,"cm"))
p2 <- p_list[[2]] + 
        ylab("")+
        xlab("-log10(q-value)")+
        ggtitle("HSC2") +
        theme(text = element_text(color="black",family="sans"),
            axis.text.x = element_text(size=8))+
        force_panelsizes(row=unit(5,"cm"),col=unit(8,"cm"))
pdf("plot/z14.hsc.marker-grp-GO.barplot.pdf")
print(plot_grid(plotlist = c(list(p1), list(p2)), ncol = 1))
dev.off()



############ GSEA
load("/home/mayao/genome_file/GeneSets/2017_NAR.2020_CSC.2021_Blood.gene_set_list.rdata")
load("/home/mayao/genome_file/GeneSets/df_geneset_ls.gs_ls.rdata")
gene_set_list_sub <- c(gene_set_list, gs_ls)

processed_gsls <- list()
for (term in names(gene_set_list_sub)){
    df <- data.frame(gs_name = term, gene_symbol=gene_set_list_sub[[term]])
    processed_gsls[[term]] <- df
}
processed_gsls_combine <- Reduce(rbind,processed_gsls)
head(processed_gsls_combine)

markers <- FindMarkers(obj, ident.1 = "HSC2", ident.2 = "HSC1")
markers <- markers[order(-markers$avg_log2FC),]
genelist <- structure(markers$avg_log2FC, names=rownames(markers))

res <- tryCatch({clusterProfiler::GSEA(genelist, TERM2GENE = processed_gsls_combine, pvalueCutoff = 1, seed = 1234)},
               error = function(e) {NULL})

res_df <- res@result[,c("ID", "NES", "p.adjust", "qvalue")]
res_df <- res_df[order(-res_df$NES),]
res_df

## *********** Figure 2D ***********
selected_terms <- c("LT_HSC", "MolO_Wilson_et_al", "Serial_engraftment_Rodriguez_et_al","ST_HSC","MPP5","MPP4")

for (id in selected_terms){
    p <- gseaNb(object = res, lineSize=1, geneSetID = id,
            segCol = "red", curveCol = rep("green",3), 
            addPval = TRUE, pvalSize = 4, pCol = "black",
            pvalX = 0.85, pvalY = 0.55
            )
    p[[1]] <- p[[1]] + 
                theme(axis.text = element_text(family = "sans", size = 12)) +
                ggtitle(id)
    ggsave(p, filename = paste0("plot/05.hsc.GseaVis.",id,".pdf"))
}






## ***************************************************************************************************************
## Combined with Previous data Sommerkamp et al.
## ***************************************************************************************************************
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

## set parameters after testing
nfeatures <- 2000
ndim <- 6
integrate_method <- "cca"
wt <- 25

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

## *********** Figure 2E, S2A ***********
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
p2 <- ggplot(umap, aes(x=!!ensym(reduction_1), y=!!ensym(reduction_2), color=!!ensym(group.by.new)))+
        geom_point(size = 0.1, alpha=0.9) +
        scale_color_manual(values = colCls, name="") +
        plotTheme +
        labs(x="UMAP 1", y="UMAP 2") +
        guides(color = guide_legend(override.aes = list(size=3))) +
        ggtitle("") +
        force_panelsizes(rows=unit(8,"cm"), cols=unit(9,"cm"))
pdf("plot/27.hsc-combine.umap-featureplot.pdf", width = 30, height = 10)
print(plot_grid(plotlist=c(list(p1), list(p2)), ncol = 3))
dev.off()

## *********** Figure 2F ***********
library(gtools)
library(ggalluvial)
hsc <- combine %>% subset(subset = Name == "hsc")
orig_df <- hsc@meta.data %>% as.data.frame()
table(orig_df$Group, orig_df$Cluster_new)

df_count <- plyr::ddply(orig_df, .(Group, Cluster_new), summarize, Counts = length(Cluster_new)) 
df_perc <- plyr::ddply(df_count, .(Group), summarize,  Cluster_new =  Cluster_new, Percent = Counts/sum(Counts)*100) 
df_perc$Cluster_new <- factor(df_count$Cluster_new, levels = rev(levels(df_count$Cluster_new)))
table(df_perc$Cluster_new)

p4 <- ggplot(df_perc, aes(x = Group, y = Percent, fill = Cluster_new, stratum = Cluster_new, alluvium =  Cluster_new)) + 
        geom_flow(width = 0.5, alpha = 0.3, knot.pos=0, color = 'white') + 
        geom_col(width = 0.5, color = 'white', alpha = 0.9) + 
        scale_y_continuous(expand = c(0, 0)) +  
        scale_fill_manual(values = rev(colCls)) +
        xlab("") + 
        ylab("Explained variation (%)") + 
        theme_classic() +
        theme(legend.title = element_blank())+
        force_panelsizes(rows=unit(8,"cm"), cols=unit(8,"cm"))
p4
ggsave("plot/28.hsc-combine.component.alluvial.pdf")


## *********** Figure S2A ***********
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
ggsave("plot/27.hsc-combine.p_module.cluster.featureplot.pdf")


## *********** Figure S2B ***********
### Markers Plots
load("Rdata/03.hsc-sk.combine.rdata")
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
ggsave(p, width = 10, height = 10, filename = "plot/30.hsc-combine.lineage-markers.dotplot.pdf")




## ***************************************************************************************************************
## HSC2 subsets
## ***************************************************************************************************************
load("Rdata/02.hsc.seuratobj.after_qc_normed.rdata")
obj <- hsc
obj <- hsc %>% subset(subset = Group == "HSC2")

## calculate mito
obj <- PercentageFeatureSet(obj, 
			pattern = "^mt-", 
			col.name = "pct_counts_Mito")

## logNorm
obj <- NormalizeData(obj, assay="RNA")

## cc
s.genes <- cc.genes$s.genes %>% str_to_title()
g2m.genes <- cc.genes$g2m.genes %>% str_to_title()
obj <- CellCycleScoring(obj , g2m.features=g2m.genes, s.features=s.genes)
obj$cc_difference <- obj$S.Score - obj$G2M.Score
obj$Phase <- factor(obj$Phase, levels=c("G1","S","G2M"))

## set parameters after testing
nfeatures <- 500
ndim <- 5
optRes <- 0.4

## Variable features 
obj <- FindVariableFeatures(obj, assay="RNA", nfeatures=nfeatures)

## scale
obj <- ScaleData(obj, assay="RNA", vars.to.regress=c('pct_counts_Mito',"cc_difference"))

## PCA
obj <- RunPCA(obj, assay="RNA", npcs=30)

## Clustering
obj <- FindNeighbors(obj, reduction = "pca" , dims = 1:ndim, verbose = FALSE)
obj <- FindClusters(obj,resolution = optRes, algorithm=4,verbose = FALSE)
table(obj$seurat_clusters)

## Rename
obj$Cluster_new <- paste0("S",obj$seurat_clusters)
obj$Cluster_new <- factor(obj$Cluster_new)
table(obj$Cluster_new)

## UMAP
reduction_name <- paste0("umap")
obj <- RunUMAP(obj, 
            reduction = "pca", 
            reduction.name = reduction_name, 
            dims=1:ndim, 
            verbose=FALSE)
## save
hsc2 <- obj
save(hsc2, file="Rdata/04.hsc2.subsets.rdata")

## Plot
load("Rdata/04.hsc2.subsets.rdata")
obj <- hsc2
## 1. UMAP showing srt clusters 
## *********** Figure 2G ***********
reduction_name <- "umap"
reduction_1 <- paste0(gsub(".", "", reduction_name, fixed=T),"_1")
reduction_2 <- paste0(gsub(".", "", reduction_name, fixed=T),"_2")
group.by.new <- "seurat_clusters"
umap <- FetchData(obj, vars = c(reduction_1, reduction_2, group.by.new))

colCls <- c("#85B3E5", "#F7AE24") 
p1 <- ggplot(umap, aes(x=!!ensym(reduction_1), y=!!ensym(reduction_2), color=!!ensym(group.by.new)))+
        geom_point(size = 0.5) +
        scale_color_manual(values = colCls, name="") +
        theme_classic() +
        labs(x="UMAP 1", y="UMAP 2") +
        guides(color = guide_legend(override.aes = list(size=3))) +
        ggtitle("") +
        theme(text = element_text(color="black")) +
        force_panelsizes(rows=unit(3,"cm"), cols=unit(3,"cm"))
p1
ggsave("plot/02.hsc-subsets.srt_cls.umap.pdf")


### 2 box plot with statistical analysis
## *********** Figure 2H ***********
load("Rdata/04.hsc2.subsets.rdata")
load("Rdata/02.hsc.seuratobj.after_qc_normed.rdata")
hsc$subsets <- "HSC1"
hsc2_cells <- rownames(hsc@meta.data[hsc$Group == "HSC2",])
hsc@meta.data[hsc2_cells,]$subsets <- hsc2@meta.data[hsc2_cells,]$Cluster_new
obj <- hsc

library(ggpubr)
library(ggbeeswarm)
library(RColorBrewer)
my_comparisons <- list(c("S1", "S2"), c("HSC1","S1"), c("HSC1", "S2"))

features <- c("Flt3", "Cd34")
boxdata <- FetchData(obj, vars=c(features,"subsets"))
boxdata$subsets <- gsub("1", "S1", boxdata$subsets)
boxdata$subsets <- gsub("2", "S2", boxdata$subsets)
boxdata$subsets <- gsub("HSCS1", "HSC1", boxdata$subsets)
boxdata$subsets <- factor(boxdata$subsets, levels=c("HSC1", "S1", "S2"))
table(boxdata$subsets)

p_val_ls <- list()
for (feature in features){
    for (n in 1:length(my_comparisons)){
        compare <- my_comparisons[[n]]
        compare_name <- paste0(compare[1], "-", compare[2])
        message("********** ", compare[1], " vs ", compare[2], " **********")
        boxdata_sub <- boxdata %>% subset(subsets %in% c(compare[1], compare[2]))
        tmp <- wilcox.test(boxdata_sub[,feature] ~ boxdata_sub[,"subsets"])
        p_val_ls[[paste0(feature, "_" , compare_name)]] <- data.frame(Feature=feature, Comparison = compare_name, Pvalue = tmp$p.value)
    }
}
p_val_df <- Reduce(rbind, p_val_ls)
p_val_df$BH_padj <- p.adjust(p_val_df$Pvalue, method = "BH")
p_val_df$Stars <- cut(p_val_df$BH_padj,
                    breaks = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                    labels = c("****", "***", "**", "*", "ns"))

p_ls <- list()
for (feature in features){
    p_val_df_sub <- p_val_df %>% subset(Feature == feature)

    y_position <- seq(from = max(boxdata[,feature])+0.2, by = 0.5, length.out = length(my_comparisons))
    y_limit <- c(0, max(y_position)+0.5)
    p_ls[[feature]] <- ggplot(boxdata, aes(x=subsets, y=!!ensym(feature), fill=subsets))+
                            geom_boxplot(color = "black", outlier.shape = NA) + 
                            scale_fill_manual(values = c("white", colCls))+
                            ggsignif::geom_signif(comparisons = my_comparisons,
                                                annotations = p_val_df_sub$Stars,  
                                                map_signif_level = TRUE, 
                                                y_position = y_position,  
                                                tip_length = 0.03, 
                                                size = 0.3, 
                                                textsize = 3, 
                                                color = "black") +
                            scale_y_continuous(limits = y_limit)+
                            labs(y="Expression",x="")+ 
                            scale_color_manual(values = c("white", colCls))+
                            theme_bw() + 
                            theme(panel.border = element_blank(),
                                panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(),
                                axis.line = element_line(colour = "black"),
                                legend.position="none", 
                                axis.text = element_text(color = "black", size = 12))+
                            ggtitle(feature)+
                            force_panelsizes(rows=unit(6,"cm"), cols=unit(5,"cm"))
}
pdf(paste0("plot/03.hsc-subsets.HSC1-added.srt_cls.CD34-Flt3.boxplot.pdf"))
print(plot_grid(plotlist = p_ls, ncol = 2))
dev.off()