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
### Define color palettes
colGroup <- c("HSC1"="black", "HSC2"="red", 
              "EMPP6"="#FF80C0","EMPP5"="#188386", "Ly-I"="#808080"
              )

## set python options and load modules
python_path <- '/home/mayao/.conda/envs/seurat/bin/'
use_python(python_path)
### set path
setwd("/home/mayao/LABdata/MY/hsc1_hsc2/mpp")
### set Seurat
options(Seurat.object.assay.version = "v5")


### ******************************************************************************************
### Pre-process
### ******************************************************************************************
meta <- openxlsx::read.xlsx("meta/I05.MPP.xlsx", sheet = 1, rowNames = FALSE)
rownames(meta) <- meta$Cell

## remove non-coding and duplicated genes
origData <- read.table("count/MPP.rawCount.txt", header = T) # 53379 291
cleanData <- keepCodingGenes(origData) %>% removeDuplicatedGenes()  # 21940 289
## only coding genes and no duplication
cleanMeta <- meta[colnames(cleanData[1:ncol(cleanData)-1]) ,] 
cleanMeta <- cleanMeta[!is.na(cleanMeta$Group),]
cleanData <- cleanData[, c(rownames(cleanMeta),"ENSEMBL")]

## save cleanData and cleanMeta
write.table(cleanMeta, file = "meta/EMPP6EMPP5LyI_sampleInfo.txt", quote = FALSE, sep = "\t", row.names = FALSE)
save(cleanData, cleanMeta, file = "Rdata/01.mpp.cleanData.cleanMeta.rdata" )


### ******************************************************************************************
### Seurat object
### ******************************************************************************************
load("Rdata/01.mpp.cleanData.cleanMeta.rdata")
dim(cleanData)
mpp <- CreateSeuratObject(counts = cleanData[,1:ncol(cleanData)-1], 
                                  meta.data = cleanMeta, 
                                  assay = "RNA", 
                                  min.cells = 3, 
                                  min.features = 200)
mpp <- PercentageFeatureSet(mpp, 
                           pattern = "^mt-", 
                           col.name = "pct_counts_Mito")
is.na(mpp$CD135_median) %>% table()

save(mpp,file = 'Rdata/temp01.mpp.seuratObj.rdata')



### ******************************************************************************************
### QC
### ******************************************************************************************
load("Rdata/temp01.mpp.seuratObj.rdata")
obj <- mpp

low_dg <- 1000 
low_umi <- 1e+04
high_dg <- 5000 
high_umi <- 5e+05
mito <- 3

obj <- subset(obj, subset=nCount_RNA > low_umi & 
                          nCount_RNA < high_umi &
                          nFeature_RNA > low_dg &
                          nFeature_RNA < high_dg &
                          pct_counts_Mito < mito)
obj$is_facs_null <- FALSE
obj@meta.data[is.na(obj$CD135_median), "is_facs_null"] <- TRUE 



### ******************************************************************************************
### Normalize and CC
### ******************************************************************************************
s.genes <- cc.genes$s.genes %>% str_to_title()
g2m.genes <- cc.genes$g2m.genes %>% str_to_title()

obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
obj <- NormalizeData(obj, assay="RNA")
obj <- FindVariableFeatures(obj, assay="RNA")

## cc
obj <- CellCycleScoring(obj, g2m.features=g2m.genes, s.features=s.genes)
obj$cc_difference <- obj$S.Score - obj$G2M.Score
obj$Phase <- factor(obj$Phase, levels=c("G1","S","G2M"))
table(obj$Phase, obj$Group)

mpp <- obj
save(mpp, file="Rdata/temp02.mpp.seurat.after_qc_norm.lognorm.rdata")



### ******************************************************************************************
### Demension reduction
### ******************************************************************************************
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

## UMAP
reduction_name <- paste0("umap")
obj <- RunUMAP(obj, 
            reduction = method, 
            reduction.name = reduction_name, 
            dims=1:ndim, 
            verbose=FALSE)

mpp <- obj
save(mpp, file = "Rdata/02.mpp.seuratobj.rdata")

## ********* Figure S5A *********
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



### ******************************************************************************************
### Markers between any 2 of the 3 groups
### ******************************************************************************************
Idents(obj) <- obj$Group
## create excel
xlsx_file <- "03.mpp.markers-between-grp.xlsx"
wb <- my_xlsx_create(xlsx_file)
## prerequisite
num <- length(levels(obj$Group))
idx <- 1
p_ls <- list()
while ( idx < num ){
    ident.1 <- levels(obj$Group)[idx]
    for (idx.2 in seq(idx, num, 1)[-1]){
        ident.2 <- levels(obj$Group)[idx.2]

        ## change ident name
        message(ident.1, " v.s. ", ident.2)

        ## markers
        res_ls <- my_markers(obj, find.all=FALSE, ident.1=ident.1, ident.2=ident.2, top_n=20, plot.title=paste0("Markers between groups\n",ident.1," vs. ", ident.2))
        markers <- res_ls[["markers_clean"]]
        p_ls[[paste0(ident.1," vs. ", ident.2)]] <- res_ls[["p"]]
        go <- res_ls[["enrich_res"]][["go"]]
        kegg <- res_ls[["enrich_res"]][["kegg"]]
        
        ## save
        wb <- my_xlsx_write(wb,sheetName=paste0(ident.1,"_", ident.2), data=markers, rowNames=TRUE)
        wb <- my_xlsx_write(wb,sheetName=paste0(ident.1,"_", ident.2,"-GO"), data=go, rowNames=TRUE)
        wb <- my_xlsx_write(wb,sheetName=paste0(ident.1,"_", ident.2,"-KEGG"), data=kegg, rowNames=TRUE)
    }
    idx <- idx + 1
    my_xlsx_save(wb, xlsx_file)
}


### ******************************************************************************************
### GSEA between any 2 of 3 groups
### ******************************************************************************************
load("/home/mayao/genome_file/GeneSets/Mm.msigdb_ls.rdata")
names(msigdb_ls)
gene_set_list <- msigdb_ls[["mh.all"]]

processed_gsls <- list()
for (term in names(gene_set_list)){
    df <- data.frame(gs_name = term, gene_symbol=gene_set_list[[term]])
    processed_gsls[[term]] <- df
}
processed_gsls_combine <- Reduce(rbind,processed_gsls)
head(processed_gsls_combine)

geneset_cat <- "Hallmarks"
Idents(obj) <- "Group"
num <- length(levels(obj$Group))
idx <- 1
res_ls <- list()
while ( idx < num ){
    ident.1 <- levels(obj$Group)[idx]
    for (idx.2 in seq(idx, num, 1)[-1]){
        ident.2 <- levels(obj$Group)[idx.2]

        markers <- FindMarkers(obj, ident.1 = ident.1, ident.2 = ident.2, logfc.threshold = 0, min.pct = 0)
        markers <- markers[order(-markers$avg_log2FC),]
        genelist <- structure(markers$avg_log2FC, names=rownames(markers))

        ## GSEA
        res <- tryCatch({GSEA(genelist, TERM2GENE = processed_gsls_combine, pvalueCutoff = 1, seed = 1234)},
                    error = function(e) {NULL})

        res_df <- res@result[,c("ID", "NES", "p.adjust", "qvalue")]
        res_df <- res_df[order(-res_df$NES),]

        ## res_ls
        res_ls[[paste0(ident.1,"-",ident.2)]] <- temp@result[,c("ID", "NES", "p.adjust", "qvalue")]
        res_ls[[paste0(ident.1,"-",ident.2)]]$Compare <- paste0(ident.1,"-",ident.2)

        ## create excel
        xlsx_file <- paste0("14.mpp.",ident.1,"-vs-",ident.2,".GSEA-Hallmarks.xlsx")
        rds_file <- paste0("Rdata/Rds06.mpp.res.",ident.1,"-vs-",ident.2,".GSEA-Hallmarks.rds")

        wb <- my_xlsx_create(xlsx_file)
        wb <- my_xlsx_write(wb,sheetName=paste0(ident.1,"-vs-",ident.2), data=res_df, rowNames=TRUE)
        my_xlsx_save(wb, xlsx_file)

        ## Save RDS 
        saveRDS(res, file = rds_file)
    }
    idx <- idx + 1
}

## ********* Figure S5B *********
res_ls_df <- Reduce(rbind, res_ls)
res_ls_df$Stars <- cut(res_ls_df$p.adjust,
                           breaks = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                           labels = c("****", "***", "**", "*", ""))
kept_term <- res_ls_df[res_ls_df$Stars != "","ID"]

ggplot(res_ls_df[res_ls_df$ID %in% kept_term, ], aes(x = Compare, y = ID, fill = NES)) +
    geom_tile(na.rm = TRUE, alpha = 0.7) +
    scale_fill_gradientn(colors = c("#0000ff", "#ffffff",  "#ff0000")) +
    geom_text(aes(label=Stars), size = 3)+
    labs(title = paste0(""), x = "", y = "") +
    theme(panel.grid = element_blank(), panel.background = element_blank(), axis.ticks = element_blank())+
    theme(axis.text = element_text(color = "black"), 
          axis.text.x = element_text(angle = 45, hjust = 1))+
    ggtitle("")+
    force_panelsizes(rows=unit(6,"cm"), cols=unit(3,"cm"))
ggsave("plot/71.mpp.GSEA-Hallmarks.tileplot.EMPP6-vs-allothers.pdf", height = 10, width = 10)




### ******************************************************************************************
### Redar plot
### ******************************************************************************************
load("/home/mayao/genome_file/GeneSets/Konturek.gs_ls.rdata")
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

## ********* Figure S5C *********
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
ggsave("plot/43.mpp.Radar.Konturek.pdf")




### ******************************************************************************************
### Combine with HSC1 and HSC2
### ******************************************************************************************
## Loading hsc data
load("/home/mayao/LABdata/MY/hsc1_hsc2/new/Rdata/02.hsc.seuratobj.after_qc_normed.rdata")
hsc$Name <- "hsc"
colnames(hsc) <- paste0("HSC.", colnames(hsc))

load("Rdata/02.mpp.seuratobj.rdata")
mpp$Name <- "mpp"
colnames(mpp) <- paste0("MPP.", colnames(mpp))

## Create combined Seurat object
counts_ls <- c(list(hsc[['RNA']]$counts), list(mpp[['RNA']]$counts))
names(counts_ls) <- c("hsc", "mpp")
meta_ls <- c(list(as.data.frame(hsc@meta.data[,c("Group",  "Name")])),
             list(as.data.frame(mpp@meta.data[,c("Group", "Name")])) )
meta <- Reduce(rbind, meta_ls)
combine <- CreateSeuratObject(counts=counts_ls,meta.data=meta)

## correct names
combine$Group <- factor(combine$Group, levels = c("HSC1", "HSC2", "EMPP6", "EMPP5", "Ly-I"))

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
                        verbose=F,
                        dims = 1:ndim,
                        seed.use = 1234
                        )

## UMAP
combine <- RunUMAP(combine, 
    reduction = "harmony", 
    reduction.name = "umap.harmony", 
    dims=1:ndim, 
    n.neighbors = 60,
    min.dist = 0.2,
    spread = 0.5,
    metric = "euclidean",
    seed.use = 1234,
    verbose=FALSE)
    
## Join the layers
combine[["RNA"]] <- JoinLayers(combine[["RNA"]])

save(combine, file = "Rdata/RCombineHSC1.02.hsc-EMPP6EMPP5LyI.combine-obj.rdata")




### ******************************************************************************************
### Trajectories
### ******************************************************************************************
load("Rdata/temp13.blood_new.before_any_process.notfiltergenes.rdata")
obj <- blood_new %>% subset(Group %in% c("HSC1", "MPP6", "MPP5", "MPP7", "MPP1", "MPP2", "MPP3","MPP4",
                                         "CMP", "CLP", "GMP", "MEP"))
obj$Group <- gsub("MPP6", "EMPP6", obj$Group)
obj$Group <- gsub("MPP5", "EMPP5", obj$Group)
obj$Group <- gsub("MPP7", "Ly-I", obj$Group)
obj$Group <- factor(obj$Group, levels = c("HSC1", "EMPP6", "EMPP5", "Ly-I", "MPP1", "MPP2", "MPP3","MPP4",
                                         "CMP", "CLP", "GMP", "MEP") )
table(obj$Group)

# pct_counts_Mito
obj <- PercentageFeatureSet(obj, 
			pattern = "^mt-", 
			col.name = "pct_counts_Mito")

## log norm
obj <- NormalizeData(obj, assay="RNA")

## cc
s.genes <- cc.genes$s.genes %>% str_to_title()
g2m.genes <- cc.genes$g2m.genes %>% str_to_title()
obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
obj <- CellCycleScoring(obj, g2m.features=g2m.genes, s.features=s.genes)
obj$cc_difference <- obj$S.Score - obj$G2M.Score
obj$Phase <- factor(obj$Phase, levels=c("G1","S","G2M"))
table(obj$Phase, obj$Group)


## Set parameters after testing
nfeatures <- 2000
ndim <- 5

## Split the layers
obj[["RNA"]] <- split(obj[["RNA"]], f=obj$Name)

## Variable features 
obj <- FindVariableFeatures(obj, assay="RNA", nfeatures=nfeatures)

## Scale
obj <- ScaleData(obj, assay="RNA", vars.to.regress=c('pct_counts_Mito',"cc_difference"))

## PCA
obj <- RunPCA(obj, assay="RNA", npcs=30)

## Integragate
set.seed(1234)
tmp_obj <- my_integrate(obj, assay="RNA", method="jpca", wt = 15, normalization.method = "LogNormalize", ndim=ndim)

## UMAP
for (pre_reduc in c("pca", "jpca")){
      reduction_name <- paste0("umap.", pre_reduc)
      tmp_obj <- RunUMAP(tmp_obj, 
            reduction = pre_reduc, 
            reduction.name = reduction_name, 
            dims=1:ndim, 
            n.components = 3,
            verbose=FALSE)
}
DimPlot(tmp_obj, group.by = "Group", reduction = "umap.jpca", pt.size = 1.5) +
      scale_color_manual(values = colGroup)

## Join the layers
tmp_obj[["RNA"]] <- JoinLayers(tmp_obj[["RNA"]])

mu_hspc_12 <- tmp_obj
saveRDS(mu_hspc_12, file = "Rdata/21.mu_hspc_12.rds")

## Diffusion map
obj <- readRDS("Rdata/21.mu_hspc_12.rds")
obs <- c("Group", "Phase")
group.by <- "Group"
pca_name <- "jpca"

return_list <- setup_anndata(obj, obs=obs, pca_name=pca_name, umap_name=paste0('umap.', pca_name), n_pcs=6, n_neighbors=30)
sc <- return_list[['sc']]
adata <- return_list[['adata']]

sc$tl$diffmap(adata)
oupDR <- py_to_r(adata$obsm['X_diffmap'])
rownames(oupDR) <- colnames(obj)
colnames(oupDR) <- paste0("DC_", 0:(15-1))
oupDR <- oupDR[, paste0("DC_", seq(6))]
obj[["diffmap"]] = CreateDimReducObject(embeddings = oupDR, key = "DC_",
                                        assay = DefaultAssay(obj))

library(plotly)
dc <- Embeddings(obj, reduction = "diffmap") %>% as.data.frame()
dc$Group <- obj$Group

p <- plot_ly(x=as.numeric(dc$DC_1), y=dc$DC_2, z=dc$DC_3, type="scatter3d", mode="markers", color=dc$Group, colors = colGroup)
htmlwidgets::saveWidget(p, file = paste0("plot/125.mu_hspc_12.3D-diffusionmap.html"))



### ******************************************************************************************
### Markers between HSC1 and EMPP6
### ******************************************************************************************
load("Rdata/RCombineHSC1.02.hsc-EMPP6EMPP5LyI.combine-obj.rdata")
obj <- combine
group.by <- "Group"
Idents(obj) <- group.by 

res_ls <- my_markers(obj, find.all=FALSE, ident.1 = "EMPP6", ident.2 = "HSC1", group.by = group.by, whether_plot = FALSE)
markers <- res_ls[["markers"]]

## box plot
## ********* Figure 5A *********
library(ggpubr)
library(ggbeeswarm)
library(RColorBrewer)
sub_obj <- obj %>% subset(Group %in% c("HSC1", "EMPP6"))

features_E6 <- c("Myl10",  "Plac8",  "Ltb",  "Tmsb10",  "Ier2",  "Ramp1",  "Cd47",  "Arpc1b",  "Ifitm1",  "Actb",  "Tmsb4x",  "Ifitm2",  "Cfl1",  "Ptma")
features_hsc1 <- c("Serinc3",  "Npm1",  "Sord",  "MmrE6",  "Acadl",  "Upp1",  "Itgb3",  "Tgm2",  "Apoe",  "Trim47",  "Itga2b",  "Clec1a",  "Aldh1a1",  "Sult1a1")

features <- c(features_E6, rev(features_hsc1))
boxdata <- FetchData(sub_obj, vars=c(features,"Group"))
boxdata_mut <- tidyr::pivot_longer(boxdata, cols = -Group, names_to = "Feature", values_to = "Expression")
boxdata_mut$Feature <- factor(boxdata_mut$Feature, levels = features)
table(boxdata_mut$Feature)

my_comparisons <- list(c("HSC1", "EMPP6"))
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
ggsave(paste0("plot/87.hsc-EMPP6EMPP5LyI.markers.EMPP6-HSC1.boxplot.pdf"), width = 15, height = 15)



### ******************************************************************************************
### GSEA
### ******************************************************************************************
load("/home/mayao/genome_file/GeneSets/KEGG.kegg_gs.mapping.pathways.rdata")
gene_set_list <- kegg_gs

gene_set_list_sub <- gene_set_list
processed_gsls <- list()
for (term in names(gene_set_list_sub)){
    df <- data.frame(gs_name = term, gene_symbol=gene_set_list_sub[[term]])
    processed_gsls[[term]] <- df
}
processed_gsls_combine <- Reduce(rbind,processed_gsls)
head(processed_gsls_combine)

load("Rdata/RCombineHSC1.02.hsc-EMPP6EMPP5LyI.combine-obj.rdata")
obj <- combine

geneset_cat <- "KEGG" 
Idents(obj) <- "Group"
ident.1 <- "EMPP6"
ident.2 <- "HSC1"

markers <- FindMarkers(obj, ident.1 = ident.1, ident.2 = ident.2, logfc.threshold = 0, min.pct = 0)
markers <- markers[order(-markers$avg_log2FC),]
genelist <- structure(markers$avg_log2FC, names=rownames(markers))

res <- tryCatch({GSEA(genelist, TERM2GENE = processed_gsls_combine, pvalueCutoff = 1, seed = 1234)},
              error = function(e) {NULL})

res_df <- res@result[,c("ID", "NES", "p.adjust", "qvalue")]
res_df <- res_df[order(-res_df$NES),]
res_df$Group <- ifelse(res_df$NES > 0, "EMPP6", "HSC1")
res_df$Group <- factor(res_df$Group, levels = c("EMPP6", "HSC1"))
res_df <- res_df[order(res_df$NES,decreasing = TRUE), ]


## Stars
res_df$Stars <- cut(res_df$p.adjust,
                           breaks = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                           labels = c("****", "***", "**", "*", ""))
kept_term <- res_df[res_df$Stars != "","ID"]
res_df <- res_df[res_df$ID %in% kept_term,]

selected_terms_EMPP6 <- c("Cytoskeleton in muscle cells ", "Rap1 signaling pathway ", "IL-17 signaling pathway ", "TNF signaling pathway ",
                    "Chemokine signaling pathway ", "NOD-like receptor signaling pathway ", "T cell receptor signaling pathway ", "Efferocytosis ", 
                    "B cell receptor signaling pathway ", "JAK-STAT signaling pathway ")
selected_terms_HSC1 <- c("Aminoacyl-tRNA biosynthesis ", "Proteasome ", "Ribosome ", "RNA polymerase ", "Nucleotide excision repair ",
                         "Mismatch repair ", "Glutathione metabolism ", "Cobalamin transport and metabolism ", "Glycerophospholipid metabolism ",
                         "Pyruvate metabolism ")

res_df_sub <- res_df[res_df$ID %in% c(selected_terms_EMPP6, selected_terms_HSC1),]
res_df_sub <- res_df_sub[order(res_df_sub$NES, decreasing = TRUE), ]
res_df_sub$ID <- factor(res_df_sub$ID, levels = rev(res_df_sub$ID))

## ********* Figure 5B *********
p1<- ggplot(res_df_sub[1:10,], aes(x = NES, y = ID, fill = Group))+
        geom_bar(stat = "identity", width = 0.5)+
        scale_fill_manual(values = colGroup) +
        geom_text(aes(label=Stars), size = 3)+
        labs(title = paste0(""), x = "", y = "") +
        theme(panel.grid = element_blank(), panel.background = element_blank())+
        theme(axis.text = element_text(color = "black"))+
        ggtitle("")+
        force_panelsizes(rows=unit(4,"cm"), cols=unit(2,"cm"))

p2 <- ggplot(res_df_sub[11:20,], aes(x = NES, y = ID, fill = Group))+
        geom_bar(stat = "identity", width = 0.5)+
        scale_fill_manual(values = colGroup) +
        geom_text(aes(label=Stars), size = 3)+
        labs(title = paste0(""), x = "", y = "") +
        theme(panel.grid = element_blank(), panel.background = element_blank())+
        theme(axis.text = element_text(color = "black"))+
        ggtitle("")+
        force_panelsizes(rows=unit(4,"cm"), cols=unit(2,"cm"))
p_ls <- c(list(p1), list(p2))
pdf("plot/78.hspc.KEGG.EMPP6-VS-HSC1.barplot.selected.pdf")
print(plot_grid(plotlist = p_ls, ncol = 2))
dev.off()





### ******************************************************************************************
### Enrichment score specificalling in hematopoietic settings
### ******************************************************************************************
load("Rdata/RCombineHSC1.02.hsc-EMPP6EMPP5LyI.combine-obj.rdata")
obj <- combine

## ********* Figure 5D, S5E, S5F *********
## DF HSPCs genelist
load("/home/mayao/genome_file/GeneSets/df_geneset_ls.gs_ls.newTermsAdded.rdata")
geneset_name <- names(gs_ls)

temp_obj <- my_AddModuleScore(obj = obj, gene_set_list = gs_ls, gsids = geneset_name)
colnames(temp_obj@meta.data)
score_df <- temp_obj@meta.data[, c(4,11:ncol(temp_obj@meta.data))] %>% as.data.frame()
score_longer <- tidyr::pivot_longer(score_df, cols = -Group, names_to = "Feature", values_to = "Score")
score_ls <- split(score_longer, f = score_longer$Feature)

## calculate p value
pval_ls <- list()
fc_ls <- list()
for (gs_id in geneset_name){
    temp_df <- score_ls[[gs_id]]

    gr1 <- temp_df[temp_df$Group == "HSC1", ]$Score
    gr2 <- temp_df[temp_df$Group == "EMPP6", ]$Score
    gr3 <- temp_df[temp_df$Group == "EMPP5", ]$Score
    gr4 <- temp_df[temp_df$Group == "Ly-I", ]$Score

    temp_vec_p <- c("EMPP6vsHSC1" = wilcox.test(gr2 , gr1, alternative = "two.sided")$p.value,
                    "EMPP5vsHSC1" = wilcox.test(gr3 , gr1, alternative = "two.sided")$p.value,
                    "LyIvsHSC1" = wilcox.test(gr4 , gr1, alternative = "two.sided")$p.value,
                    "EMPP5vsEMPP6" = wilcox.test(gr3 , gr2, alternative = "two.sided")$p.value,
                    "LyIvsEMPP6" = wilcox.test(gr4 , gr2, alternative = "two.sided")$p.value,
                    "LyIvsEMPP5" = wilcox.test(gr4 , gr3, alternative = "two.sided")$p.value)

    ## calculate log2FC
    temp_vec_fc <- c("EMPP6vsHSC1" = mean(gr2) -  mean(gr1),
                     "EMPP5vsHSC1" = mean(gr3) -  mean(gr1),
                     "LyIvsHSC1" = mean(gr4) -  mean(gr1),
                     "EMPP5vsEMPP6" = mean(gr3) -  mean(gr2),
                     "LyIvsEMPP6" = mean(gr4) -  mean(gr2),
                     "LyIvsEMPP5" = mean(gr4) -  mean(gr3))
    temp_vec_fc <- signif(temp_vec_fc, digits = 2)

    pval_ls[[gs_id]] <- temp_vec_p
    if ( !is.null(temp_vec_fc) ){
        fc_ls[[gs_id]] <- temp_vec_fc
    }
}
pval_df <- Reduce(rbind, pval_ls) %>% as.data.frame()
rownames(pval_df) <- names(pval_ls)
pval_df$Feature <- rownames(pval_df)
pval_df_longer <- tidyr::pivot_longer(pval_df, cols = 1:6, names_to = "Comparison", values_to = "pval")
head(pval_df_longer)

## subset: keep only HSC1 and EMPP6
pval_df_sub <- pval_df_longer %>% subset(Comparison == "EMPP6vsHSC1")

## p value correction
pval_df_sub$padj <- p.adjust(pval_df_sub$pval, method = "BH")
pval_df_sub

## Stars
pval_df_sub$stars <- cut(pval_df_sub$padj,
            breaks = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
            labels = c("****", "***", "**", "*", "ns"))
pval_df_sub

## plot boxplot: all geneset terms
p_ls <- list()
for (gs_id in pval_df_sub$Feature){
    plot_data <- score_longer %>% subset(Feature == gs_id & Group %in% c("HSC1","EMPP6"))
    comparisons <- list(c("EMPP6", "HSC1"))

    y_position <- seq(from = max(plot_data$Score)+0.02, by = 0.1, length.out = length(comparisons))
    y_limit <- c(0, max(y_position)+0.1)

    p_ls[[paste0(gs_id)]] <- ggplot(plot_data, aes(x = Group, y = Score, color = Group)) +
                    geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.8) +
                    ggbeeswarm::geom_quasirandom(aes(color=Group),size = 0.3, alpha = 0.8)+
                    ggsignif::geom_signif(comparisons = comparisons,
                                    annotations = pval_df_sub[pval_df_sub$Feature == gs_id,]$stars,  
                                    map_signif_level = TRUE, 
                                    y_position = y_position,  
                                    tip_length = 0.03, 
                                    size = 0.3, 
                                    textsize = 3, 
                                    color = "black") +
                    scale_color_manual(values = colGroup) +
                    scale_shape_manual(values = c(16, 4))+
                    labs(title = paste0(gs_id), x = "", y = "Module score") +
                    scale_y_continuous(limits = y_limit)+
                    theme_classic()+
                    theme(axis.text = element_text(color = "black", size = 14), 
                            # axis.text.x = element_text(angle = 45, hjust = 1), 
                            plot.title = element_text(color = "black", face = "italic", size = 16), 
                            axis.title = element_text(color = "black", size = 16))+
                    force_panelsizes(rows=unit(4,"cm"), cols=unit(4,"cm"))
}
pdf("plot/120.hsc-EMPP6EMPP5LyI.df_hspc_gs.addmodulescore.boxplot.wilcoxon-BH.pdf", width = 25, height = 30)
print(plot_grid(plotlist=p_ls, ncol = 6))
dev.off()

## plot boxplot: only significant geneset terms
pval_df_sub_sig <- pval_df_sub %>% subset(padj < 0.05)
p_ls <- list()
for (gs_id in pval_df_sub_sig$Feature){
    plot_data <- score_longer %>% subset(Feature == gs_id & Group %in% c("HSC1","EMPP6"))
    comparisons <- list(c("EMPP6", "HSC1"))

    y_position <- seq(from = max(plot_data$Score)+0.02, by = 0.1, length.out = length(comparisons))
    y_limit <- c(0, max(y_position)+0.1)

    p_ls[[paste0(gs_id)]] <- ggplot(plot_data, aes(x = Group, y = Score, color = Group)) +
                    geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.8) +
                    ggbeeswarm::geom_quasirandom(aes(color=Group),size = 0.3, alpha = 0.8)+
                    ggsignif::geom_signif(comparisons = comparisons,
                                    annotations = pval_df_sub_sig[pval_df_sub_sig$Feature == gs_id,]$stars,  
                                    map_signif_level = TRUE, 
                                    y_position = y_position,  
                                    tip_length = 0.03, 
                                    size = 0.3, 
                                    textsize = 3, 
                                    color = "black") +
                    scale_color_manual(values = colGroup) +
                    scale_shape_manual(values = c(16, 4))+
                    labs(title = paste0(gs_id), x = "", y = "Module score") +
                    scale_y_continuous(limits = y_limit)+
                    theme_classic()+
                    theme(axis.text = element_text(color = "black", size = 14), 
                            # axis.text.x = element_text(angle = 45, hjust = 1), 
                            plot.title = element_text(color = "black", face = "italic", size = 16), 
                            axis.title = element_text(color = "black", size = 16))+
                    force_panelsizes(rows=unit(4,"cm"), cols=unit(4,"cm"))
}
pdf("plot/121.hsc-EMPP6EMPP5LyI.df_hspc_gs.SELECTED.addmodulescore.boxplot.wilcoxon-BH.pdf", width = 25, height = 20)
print(plot_grid(plotlist=p_ls, ncol = 6))
dev.off()


## ********* Figure 5E *********
## 2025 Cell Research
load("/home/mayao/genome_file/GeneSets/2025_CR_Zhangyi.gs_ls.rdata")
gs_ls <- c(gs_ls[c(1,3)], gs_ls[[2]])
geneset_name <- names(gs_ls)

temp_obj <- my_AddModuleScore(obj = obj, gene_set_list = gs_ls, gsids = geneset_name)
colnames(temp_obj@meta.data)
score_df <- temp_obj@meta.data[, c(4,11:ncol(temp_obj@meta.data))] %>% as.data.frame()
score_longer <- tidyr::pivot_longer(score_df, cols = -Group, names_to = "Feature", values_to = "Score")
score_ls <- split(score_longer, f = score_longer$Feature)

## calculate p value
pval_ls <- list()
fc_ls <- list()
for (gs_id in geneset_name){
    temp_df <- score_ls[[gs_id]]

    gr1 <- temp_df[temp_df$Group == "HSC1", ]$Score
    gr2 <- temp_df[temp_df$Group == "EMPP6", ]$Score
    gr3 <- temp_df[temp_df$Group == "EMPP5", ]$Score
    gr4 <- temp_df[temp_df$Group == "Ly-I", ]$Score

    temp_vec_p <- c("EMPP6vsHSC1" = wilcox.test(gr2 , gr1, alternative = "two.sided")$p.value,
                    "EMPP5vsHSC1" = wilcox.test(gr3 , gr1, alternative = "two.sided")$p.value,
                    "LyIvsHSC1" = wilcox.test(gr4 , gr1, alternative = "two.sided")$p.value,
                    "EMPP5vsEMPP6" = wilcox.test(gr3 , gr2, alternative = "two.sided")$p.value,
                    "LyIvsEMPP6" = wilcox.test(gr4 , gr2, alternative = "two.sided")$p.value,
                    "LyIvsEMPP5" = wilcox.test(gr4 , gr3, alternative = "two.sided")$p.value)

    ## calculate log2FC
    temp_vec_fc <- c("EMPP6vsHSC1" = mean(gr2) -  mean(gr1),
                     "EMPP5vsHSC1" = mean(gr3) -  mean(gr1),
                     "LyIvsHSC1" = mean(gr4) -  mean(gr1),
                     "EMPP5vsEMPP6" = mean(gr3) -  mean(gr2),
                     "LyIvsEMPP6" = mean(gr4) -  mean(gr2),
                     "LyIvsEMPP5" = mean(gr4) -  mean(gr3))
    temp_vec_fc <- signif(temp_vec_fc, digits = 2)

    pval_ls[[gs_id]] <- temp_vec_p
    if ( !is.null(temp_vec_fc) ){
        fc_ls[[gs_id]] <- temp_vec_fc
    }
}
pval_df <- Reduce(rbind, pval_ls) %>% as.data.frame()
rownames(pval_df) <- names(pval_ls)
pval_df$Feature <- rownames(pval_df)
pval_df_longer <- tidyr::pivot_longer(pval_df, cols = 1:6, names_to = "Comparison", values_to = "pval")
head(pval_df_longer)

## subset: keep only HSC1 and EMPP6
pval_df_sub <- pval_df_longer %>% subset(Comparison == "EMPP6vsHSC1")

## p value correction
pval_df_sub$padj <- p.adjust(pval_df_sub$pval, method = "BH")
pval_df_sub

## Stars
pval_df_sub$stars <- cut(pval_df_sub$padj,
            breaks = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
            labels = c("****", "***", "**", "*", "ns"))
pval_df_sub

## plot boxplot: all geneset terms
p_ls <- list()
for (gs_id in pval_df_sub$Feature){
    plot_data <- score_longer %>% subset(Feature == gs_id & Group %in% c("HSC1","EMPP6"))
    comparisons <- list(c("EMPP6", "HSC1"))

    y_position <- seq(from = max(plot_data$Score)+0.02, by = 0.1, length.out = length(comparisons))
    y_limit <- c(0, max(y_position)+0.1)

    p_ls[[paste0(gs_id)]] <- ggplot(plot_data, aes(x = Group, y = Score, color = Group)) +
                    geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.8) +
                    ggbeeswarm::geom_quasirandom(aes(color=Group),size = 0.3, alpha = 0.8)+
                    ggsignif::geom_signif(comparisons = comparisons,
                                    annotations = pval_df_sub[pval_df_sub$Feature == gs_id,]$stars,  
                                    map_signif_level = TRUE, 
                                    y_position = y_position,  
                                    tip_length = 0.03, 
                                    size = 0.3, 
                                    textsize = 3, 
                                    color = "black") +
                    scale_color_manual(values = colGroup) +
                    scale_shape_manual(values = c(16, 4))+
                    labs(title = paste0(gs_id), x = "", y = "Module score") +
                    scale_y_continuous(limits = y_limit)+
                    theme_classic()+
                    theme(axis.text = element_text(color = "black", size = 14), 
                            # axis.text.x = element_text(angle = 45, hjust = 1), 
                            plot.title = element_text(color = "black", face = "italic", size = 16), 
                            axis.title = element_text(color = "black", size = 16))+
                    force_panelsizes(rows=unit(4,"cm"), cols=unit(4,"cm"))
}
pdf("plot/122.hsc-EMPP6EMPP5LyI.zhangyi.addmodulescore.boxplot.wilcoxon-BH.pdf", width = 18, height = 10)
print(plot_grid(plotlist=p_ls, ncol = 4))
dev.off()

### combine q1,q2,q3,q4
selected <- c("q1","q2","q3","q4")
plot_data <- score_longer %>% subset(Feature %in% selected & Group %in% c("HSC1","EMPP6"))
ggplot(plot_data, aes(x = Group, y = Score, color = Group)) +
    geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.8) +
    ggbeeswarm::geom_quasirandom(aes(color = Group),size = 0.3, alpha = 0.65)+
    facet_wrap(~Feature, ncol=4)+
    scale_color_manual(values = colGroup) +
    scale_shape_manual(values = c(16, 4))+
    labs(x = "", y = "Module score") +
    scale_y_continuous(limits = c(0,2))+
    theme_classic()+
    theme(axis.text = element_text(color = "black", size = 8), 
        axis.ticks = element_line(color = "black"), 
        plot.title = element_text(color = "black", face = "italic", size = 16), 
        axis.title = element_text(color = "black", size = 8))+
    force_panelsizes(rows=unit(4,"cm"), cols=unit(3.5,"cm"))
ggsave("plot/123.hsc-EMPP6EMPP5LyI.zhangyi.q1234.addmodulescore.boxplot.wilcoxon-BH.pdf")