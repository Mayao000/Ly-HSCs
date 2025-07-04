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
library(patchwork)  %>% suppressMessages()
library(ggh4x)  %>% suppressMessages()
library(Matrix) %>% suppressMessages()
library(ggrepel) %>% suppressMessages()
library(data.table) %>% suppressMessages()
library(GSVA) %>% suppressMessages()
library(tidyverse) %>% suppressMessages()
library(clusterProfiler) %>% suppressMessages()
library(ggsci) %>% suppressMessages()
library(GseaVis) %>% suppressMessages()

col_gsea1<-pal_simpsons()(16)


### set path
setwd("/home/mayao/LABdata/MY/hsc1_hsc2/new")

### set Seurat
options(Seurat.object.assay.version = "v5")


### ----------------------
### 1. Prepare gene sets
### ----------------------
# load("/home/mayao/genome_file/GeneSets/2017_NAR.2020_CSC.2021_Blood.gene_set_list.rdata")
# gene_set_list_sub <- gene_set_list

load("/home/mayao/genome_file/GeneSets/df_geneset_ls.gs_ls.rdata")
gene_set_list_sub <- gs_ls


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
load("Rdata/02.hsc.seuratobj.after_qc_normed.rdata")
obj <- hsc
markers <- FindMarkers(obj, ident.1 = "HSC2", ident.2 = "HSC1")
markers <- markers[order(-markers$avg_log2FC),]
genelist <- structure(markers$avg_log2FC, names=rownames(markers))


res <- tryCatch({GSEA(genelist, TERM2GENE = processed_gsls_combine, pvalueCutoff = 1, seed = 1234)},
               error = function(e) {NULL})

res_df <- res@result[,c("ID", "NES", "p.adjust", "qvalue")]
res_df <- res_df[order(-res_df$NES),]
res_df


## Visualise using GseaVis
# selected_terms <- c("LT_HSC","ST_HSC","MPP5","MPP4")
selected_terms <- c("MolO_Wilson_et_al", "Serial_engraftment_Rodriguez_et_al")

for (id in selected_terms){
  p <- gseaNb(object = res, lineSize=1, geneSetID = id,
          segCol = "red", curveCol = rep("green",3), 
          addPval = TRUE, pvalSize = 4, pCol = "black",
          pvalX = 0.85, pvalY = 0.55
          )
  p[[1]] <- p[[1]] + 
              theme(axis.text = element_text(family = "sans", size = 12)) +
              ggtitle(id)
#   ggsave(p, filename = paste0("plot/04.hsc.GseaVis.Pei-Sommerkamp.",id,".pdf"))
    ggsave(p, filename = paste0("plot/05.hsc.GseaVis.DF.",id,".pdf"))
}
