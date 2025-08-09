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
library(ggpubr) %>% suppressMessages()
library(ggbeeswarm) %>% suppressMessages()

### set plot options
### Define color palettes and plot themes
colLibrary <- colorRampPalette(brewer.pal(n = 5, name = "Spectral"))(5)
colGEX <- c("grey85", brewer.pal(7, "Reds"))      # color for gene expression
colCcy <- c("black", "blue", "darkorange")        # color for cellcycle phase
plotTheme <- theme_classic(base_size = 14)
# colGroup <- colorRampPalette(brewer.pal(n = 4, name = "Spectral"))(2)
colGroup <- c("#FF80C0", "#188386", "#808080")
show_col(colGroup)
colTissue <- brewer.pal(4,'RdGy')[c(1,4)]
colCondition <- brewer.pal(4,'PuOr')[c(3,4)]

## set other options after being decided
colCls <- colorRampPalette(brewer.pal(n = 3, name = "Set1"))(3)

## set python options and load modules
python_path <- '/home/mayao/.conda/envs/seurat/bin/'
use_python(python_path)
# py_run_string("import numpy as np")
# py_run_string("import matplotlib.pyplot as pl")



### set path
setwd("/home/mayao/LABdata/MY/hsc1_hsc2/mpp")

### set Seurat
options(Seurat.object.assay.version = "v5")




## --------------------------------------------------------------------------------------------
load("Rdata/02.mpp.seuratobj.rdata")
obj <- mpp
## Change names
obj$Group <- gsub("MPP6", "N1", obj$Group)
obj$Group <- gsub("MPP5", "N2", obj$Group)
obj$Group <- gsub("MPP7", "Ly-I", obj$Group)
obj$Group <- factor(obj$Group, levels = c('N1', "N2", 'Ly-I'))
colGroup <- setNames(colGroup, nm = c('N1', "N2", 'Ly-I' ))
table(obj$Group)


## --------------------------------------------------------------------------------------------------
## Markers between groups
Idents(obj) <- "Group"

## N1-N2, N1-N3
ident.1 <- "N1"
for (ident.2 in c("N2", "Ly-I")){
    markers <- FindMarkers(obj, ident.1=ident.1, ident.2=ident.2, logfc.threshold=0, , min.pct = 0, recorrect_umi = FALSE) %>% 
                    # subset(subset=p_val_adj < 0.05) %>% 
                    plotly::arrange(desc(avg_log2FC))
    markers$gene <- rownames(markers)

    rds_file <- paste0("Rdata/Rds01.mpp.markers.no-p-cutoff.",ident.1,"-vs-",ident.2,".rds")
    saveRDS(markers, file = rds_file)


    # clean markers
    # kept_genes <- rownames(markers)[!grepl("Rpl|Rps|Gm|mt-", rownames(markers))]
    # markers_clean <-markers[kept_genes,]

    # rds_file <- paste0("Rdata/Rds15.mpp.markers.",ident.1,"-vs-",ident.2,".rds")
    # saveRDS(markers_clean, file = rds_file)

}

## N2-N3
markers <- FindMarkers(obj, ident.1="N2", ident.2="Ly-I", logfc.threshold=0, , min.pct = 0, recorrect_umi = FALSE) %>% 
                    # subset(subset=p_val_adj < 0.05) %>% 
                    plotly::arrange(desc(avg_log2FC))
markers$gene <- rownames(markers)
rds_file <- paste0("Rdata/Rds02.mpp.markers.no-p-cutoff.N2-vs-Ly-I.rds")
saveRDS(markers, file = rds_file)

# # clean markers
# kept_genes <- rownames(markers)[!grepl("Rpl|Rps|Gm|mt-", rownames(markers))]
# markers_clean <-markers[kept_genes,]

# rds_file <- paste0("Rdata/Rds15.mpp.markers.N2-vs-N3.rds")
# saveRDS(markers_clean, file = rds_file)



## load markers
res_ls <- list()
ident.1 <- "N1"
for (ident.2 in c("N2", "Ly-I")){
    rds_file <- paste0("Rdata/Rds01.mpp.markers.no-p-cutoff.",ident.1,"-vs-",ident.2,".rds")
    res_ls[[ident.2]] <- readRDS(rds_file)
    res_ls[[ident.2]]$Compare <- paste0(ident.1,"-", ident.2)
}

res_ls[["N2-Ly-I"]] <- readRDS("Rdata/Rds02.mpp.markers.no-p-cutoff.N2-vs-Ly-I.rds")
res_ls[["N2-Ly-I"]]$Compare <- "N2-Ly-I"

res_ls_df <- Reduce(rbind, res_ls)
head(res_ls_df)
## Stars
res_ls_df$Stars <- cut(res_ls_df$p_val_adj,
                           breaks = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                           labels = c("****", "***", "**", "*", ""))
kept_gene <- res_ls_df[res_ls_df$Stars != "","gene"]
kept_gene


res_ls_df[res_ls_df$Stars != "", "Compare"] %>% table()

## boxplot consistent with HSC1 vs. N1
library(ggpubr)
library(ggbeeswarm)
library(RColorBrewer)

colGroup_plot <- setNames(colGroup, nm = c("N1", "N2", "Ly-I"))
p_ls <- list()
for (ident.1 in c("N1", "N2")){
    if (ident.1 == "N1"){
        ident2_ls <- c("N2", "Ly-I")
    }else{
        ident2_ls <- "Ly-I"
    }

    for (ident.2 in ident2_ls){
        ## features to plot
        sub_df <- res_ls_df[res_ls_df$Stars != "" & res_ls_df$Compare == paste0(ident.1,"-", ident.2), ] 
        colnames(sub_df)[which(colnames(sub_df) == "gene")] <- "Feature"
        feature_1 <- sub_df[sub_df$avg_log2FC > 0, ]$Feature
        feature_2 <- sub_df[sub_df$avg_log2FC < 0, ]$Feature
        features <- c(feature_1, feature_2)

        ## subset obj
        sub_obj <- obj %>% subset(Group %in% c(ident.1, ident.2))

        ## boxplot data
        boxdata <- FetchData(sub_obj, vars=c(features,"Group"))
        boxdata_mut <- tidyr::pivot_longer(boxdata, cols = -Group, names_to = "Feature", values_to = "Expression")
        boxdata_mut$Feature <- factor(boxdata_mut$Feature, levels = features)
        table(boxdata_mut$Feature)

        ## set comparison
        my_comparisons <- list(c(ident.1, ident.2))

        ## get p value anno data
        boxdata_anno <- merge(boxdata_mut, sub_df[, c("Feature", "Stars")], by = "Feature")
        anno_1 <- boxdata_anno %>% subset(Feature %in% feature_1)
        anno_1$y_pos <- max(anno_1$Expression) - 0.05 * max(anno_1$Expression)
        anno_2 <- boxdata_anno %>% subset(Feature %in% feature_2)
        anno_2$y_pos <- max(anno_2$Expression) - 0.05 * max(anno_2$Expression)


        ## plot
        compare_name <- paste0(ident.1, "-", ident.2) 
        p_ls[[paste0(compare_name, "-1")]] <- ggplot(boxdata_mut %>% subset(Feature %in% feature_1), aes(x=Group, y=Expression, fill=Group))+
                                                geom_boxplot(color = "black", outlier.shape = NA, alpha = 1) + 
                                                facet_wrap( vars(Feature), ncol = 10 )+
                                                scale_fill_manual(values = colGroup_plot[c(ident.1, ident.2)] )+
                                                geom_text(data = anno_1, aes(x = 1.5, y = y_pos, label = Stars), inherit.aes = FALSE, size = 6, color = "black") +
                                                labs(y="Expression",x="")+ 
                                                scale_color_manual(values = colGroup_plot[c(ident.1, ident.2)] )+
                                                theme_bw() + 
                                                theme(
                                                        # panel.border = element_blank(),
                                                        panel.grid.major = element_blank(),
                                                        panel.grid.minor = element_blank(),
                                                        axis.line = element_line(colour = "black"),
                                                        legend.position="none", 
                                                        axis.text = element_text(color = "black", size = 16),
                                                        strip.text = element_text(color = "black", size = 16, face = "italic"),
                                                        axis.title = element_text(color = "black", size = 16) )+
                                                force_panelsizes(rows=unit(3,"cm"), cols=unit(3,"cm"))
        p_ls[[paste0(compare_name, "-2")]] <- ggplot(boxdata_mut %>% subset(Feature %in% feature_2), aes(x=Group, y=Expression, fill=Group))+
                                                geom_boxplot(color = "black", outlier.shape = NA, alpha = 1) + 
                                                facet_wrap( vars(Feature), ncol = 10 )+
                                                scale_fill_manual(values = colGroup_plot[c(ident.1, ident.2)] )+
                                                geom_text(data = anno_2, aes(x = 1.5, y = y_pos, label = Stars), inherit.aes = FALSE, size = 6, color = "black") +
                                                labs(y="Expression",x="")+ 
                                                scale_color_manual(values = colGroup_plot[c(ident.1, ident.2)] )+
                                                theme_bw() + 
                                                theme(
                                                        # panel.border = element_blank(),
                                                        panel.grid.major = element_blank(),
                                                        panel.grid.minor = element_blank(),
                                                        axis.line = element_line(colour = "black"),
                                                        legend.position="none", 
                                                        axis.text = element_text(color = "black", size = 16),
                                                        strip.text = element_text(color = "black", size = 16, face = "italic"),
                                                        axis.title = element_text(color = "black", size = 16) )+
                                                force_panelsizes(rows=unit(3,"cm"), cols=unit(3,"cm"))

    }
}
pdf("plot/88.mpp.n1n2n3.revised_boxplot.markers-between-grp.pdf", height = 25, width = 15)
print(plot_grid(plotlist = p_ls, ncol = 1))
dev.off()
