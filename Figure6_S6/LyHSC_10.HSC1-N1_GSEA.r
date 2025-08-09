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





### ----------------------
### 1. Prepare gene sets
### ----------------------
# load("/home/mayao/genome_file/GeneSets/2017_NAR.2020_CSC.2021_Blood.gene_set_list.rdata")
# load("/home/mayao/genome_file/GeneSets/Reactome.reactome_gs.rdata")
# gene_set_list <- reactome_gs
# load("/home/mayao/genome_file/GeneSets/Mm.msigdb_ls.rdata")
# names(msigdb_ls)
# gene_set_list <- msigdb_ls[["m5.go.bp"]]
# gene_set_list <- msigdb_ls[["mh.all"]]

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


### ---------------------------
### 2. Prepare ranked gene list
### ---------------------------
load("Rdata/RCombineHSC1.02.hsc-n1n2n3.combine-obj.rdata")
obj <- combine

## Correct names
obj$Group <- gsub("MPP6", "N1", obj$Group)
obj$Group <- gsub("MPP5", "N2", obj$Group)
obj$Group <- gsub("MPP7", "N3", obj$Group)
obj$Group <- factor(obj$Group, levels = c("HSC1", 'N1', "N2", 'N3', paste0("MPP", seq(4))))
colGroup <- setNames(colGroup, nm = c("HSC1", paste0("MPP", seq(4)), 'N2', "N1", 'N3' )  )
table(obj$Group)


geneset_cat <- "KEGG" # Reactome GOBP Hallmarks KEGG GOBP-selected
Idents(obj) <- "Group"
ident.1 <- "N1"
ident.2 <- "HSC1"

markers <- FindMarkers(obj, ident.1 = ident.1, ident.2 = ident.2, logfc.threshold = 0, min.pct = 0)
markers <- markers[order(-markers$avg_log2FC),]
genelist <- structure(markers$avg_log2FC, names=rownames(markers))

res <- tryCatch({GSEA(genelist, TERM2GENE = processed_gsls_combine, pvalueCutoff = 1, seed = 1234)},
              error = function(e) {NULL})

res_df <- res@result[,c("ID", "NES", "p.adjust", "qvalue")]
res_df <- res_df[order(-res_df$NES),]


## create excel
if (geneset_cat == "Reactome"){
  xlsx_file <- paste0("09.hspc.",ident.1,"-vs-",ident.2,".GSEA-",geneset_cat,".xlsx")
  rds_file <- paste0("Rdata/Rds01.hspc-combine.res.",ident.1,"-vs-",ident.2,".GSEA-",geneset_cat,".rds")
}else if (geneset_cat == "GOBP"){
  xlsx_file <- paste0("10.hspc.",ident.1,"-vs-",ident.2,".GSEA-",geneset_cat,".xlsx")
  rds_file <- paste0("Rdata/Rds02.hspc-combine.res.",ident.1,"-vs-",ident.2,".GSEA-",geneset_cat,".rds")
}else if (geneset_cat == "Hallmarks"){
  xlsx_file <- paste0("11.hspc.",ident.1,"-vs-",ident.2,".GSEA-",geneset_cat,".xlsx")
  rds_file <- paste0("Rdata/Rds03.hspc-combine.res.",ident.1,"-vs-",ident.2,".GSEA-",geneset_cat,".rds")
}else if (geneset_cat == "KEGG"){
  xlsx_file <- paste0("16.hspc.",ident.1,"-vs-",ident.2,".GSEA-",geneset_cat,".xlsx")
  rds_file <- paste0("Rdata/Rds07.hspc-combine.res.",ident.1,"-vs-",ident.2,".GSEA-",geneset_cat,".rds")
}else{
  xlsx_file <- paste0("Selected.hspc.",ident.1,"-vs-",ident.2,".GSEA-",geneset_cat,".xlsx")
  rds_file <- paste0("Rdata/RdsSelected.hspc-combine.res.",ident.1,"-vs-",ident.2,".GSEA-",geneset_cat,".rds")
}
if (file.exists(paste0("table/", xlsx_file))){
  file.remove(paste0("table/", xlsx_file))
}
wb <- my_xlsx_create(xlsx_file)
wb <- my_xlsx_write(wb,sheetName=paste0(ident.1,"-vs-",ident.2), data=res_df, rowNames=TRUE)
my_xlsx_save(wb, xlsx_file)

## Save RDS 
saveRDS(res, file = rds_file)
## -----------------------------------------------------------------------------------------------------------------





## -----------------------------------------------------------------------------------------------------------------
## Visualize --- barplot
geneset_cat <- "KEGG"
ident.1 <- "N1"
ident.2 <- "HSC1"
if (geneset_cat == "Reactome"){
    rds_file <- paste0("Rdata/Rds01.hspc-combine.res.",ident.1,"-vs-",ident.2,".GSEA-Reactome.rds")
  }else if (geneset_cat == "GOBP"){
    rds_file <- paste0("Rdata/Rds02.hspc-combine.res.",ident.1,"-vs-",ident.2,".GSEA-GOBP.rds")
  }else if (geneset_cat == "Hallmarks"){
    rds_file <- paste0("Rdata/Rds03.hspc-combine.res.",ident.1,"-vs-",ident.2,".GSEA-Hallmarks.rds")
  }else if (geneset_cat == "KEGG"){
    rds_file <- paste0("Rdata/Rds07.hspc-combine.res.",ident.1,"-vs-",ident.2,".GSEA-",geneset_cat,".rds")
}else{
    rds_file <- paste0("Rdata/RdsSelected.hspc-combine.res.",ident.1,"-vs-",ident.2,".GSEA-",geneset_cat,".rds")
}
temp <- readRDS(rds_file)
res_ls_df <- temp@result[,c("ID", "NES", "p.adjust", "qvalue")]
res_ls_df$Group <- ifelse(res_ls_df$NES > 0, "N1", "HSC1")
res_ls_df$Group <- factor(res_ls_df$Group, levels = c("N1", "HSC1"))
res_ls_df <- res_ls_df[order(res_ls_df$NES,decreasing = TRUE), ]


## Stars
res_ls_df$Stars <- cut(res_ls_df$p.adjust,
                           breaks = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                           labels = c("****", "***", "**", "*", ""))
kept_term <- res_ls_df[res_ls_df$Stars != "","ID"]
res_ls_df <- res_ls_df[res_ls_df$ID %in% kept_term,]
hist(res_ls_df$NES)
table(res_ls_df$Group)

selected_terms_N1 <- c("Cytoskeleton in muscle cells ", "Rap1 signaling pathway ", "IL-17 signaling pathway ", "TNF signaling pathway ",
                    "Chemokine signaling pathway ", "NOD-like receptor signaling pathway ", "T cell receptor signaling pathway ", "Efferocytosis ", 
                    "B cell receptor signaling pathway ", "JAK-STAT signaling pathway ")
selected_terms_HSC1 <- c("Aminoacyl-tRNA biosynthesis ", "Proteasome ", "Ribosome ", "RNA polymerase ", "Nucleotide excision repair ",
                         "Mismatch repair ", "Glutathione metabolism ", "Cobalamin transport and metabolism ", "Glycerophospholipid metabolism ",
                         "Pyruvate metabolism ")

res_ls_df_sub <- res_ls_df[res_ls_df$ID %in% c(selected_terms_N1, selected_terms_HSC1),]
res_ls_df_sub <- res_ls_df_sub[order(res_ls_df_sub$NES, decreasing = TRUE), ]
res_ls_df_sub$ID <- factor(res_ls_df_sub$ID, levels = rev(res_ls_df_sub$ID))

p1<- ggplot(res_ls_df_sub[1:10,], aes(x = NES, y = ID, fill = Group))+
        geom_bar(stat = "identity", width = 0.5)+
        scale_fill_manual(values = colGroup) +
        geom_text(aes(label=Stars), size = 3)+
        labs(title = paste0(""), x = "", y = "") +
        theme(panel.grid = element_blank(), panel.background = element_blank())+
        theme(axis.text = element_text(color = "black"))+
        ggtitle("")+
        force_panelsizes(rows=unit(4,"cm"), cols=unit(2,"cm"))

p2 <- ggplot(res_ls_df_sub[11:20,], aes(x = NES, y = ID, fill = Group))+
        geom_bar(stat = "identity", width = 0.5)+
        scale_fill_manual(values = colGroup) +
        geom_text(aes(label=Stars), size = 3)+
        labs(title = paste0(""), x = "", y = "") +
        theme(panel.grid = element_blank(), panel.background = element_blank())+
        theme(axis.text = element_text(color = "black"))+
        ggtitle("")+
        force_panelsizes(rows=unit(4,"cm"), cols=unit(2,"cm"))
p_ls <- c(list(p1), list(p2))
pdf("plot/78.hspc.KEGG.N1-VS-HSC1.barplot.selected.pdf")
print(plot_grid(plotlist = p_ls, ncol = 2))
dev.off()
## -----------------------------------------------------------------------------------------------------------------








