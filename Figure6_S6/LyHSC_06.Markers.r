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
## Correct names
obj$Group <- gsub("MPP6", "N1", obj$Group)
obj$Group <- gsub("MPP5", "N2", obj$Group)
obj$Group <- gsub("MPP7", "Ly-I", obj$Group)
obj$Group <- factor(obj$Group, levels = c('N1', "N2", 'Ly-I'))
colGroup <- setNames(colGroup, nm = c('N1', "N2", 'Ly-I' ))
table(obj$Group)

## --------------------------------------------------------------------------------------------
## 1. IFN family
features <- c("Ifngr1","Ifngr2","Ifnar1","Ifnar2", "Tnfrsf1a", "Tnfrsf1b", "Fas", "Cd40", "Cd27", "Tnfrsf11a")
DotPlot(obj, features = features, group.by = "Group", scale = FALSE) +
    coord_flip()


df <- FetchData(obj, vars = c("Group", features))
df_long <- tidyr::pivot_longer(df, cols = -Group, names_to = "Feature", values_to = "Expression")
head(df_long)

ggplot(df, aes(Group, features))
## --------------------------------------------------------------------------------------------




# ## --------------------------------------------------------------------------------------------------
# ## 2. Process only HSC & MPP
# load("Rdata/temp13.blood_new.before_any_process.notfiltergenes.rdata")
# obj <- blood_new %>% subset(subset = Group %in% c("HSC1",  paste0("MPP", seq(7))))
# obj$Group <- factor(obj$Group, levels = c("HSC1", paste0("MPP", seq(7))))
# table(obj$Group)
# #  HSC1 MPP1 MPP2 MPP3 MPP4 MPP5 MPP6 MPP7 
# #  159   72   72   76   78   78   81   80


# ## pct_counts_Mito
# obj <- PercentageFeatureSet(obj, 
# 			pattern = "^mt-", 
# 			col.name = "pct_counts_Mito")

# ## log norm
# s.genes <- cc.genes$s.genes %>% str_to_title()
# g2m.genes <- cc.genes$g2m.genes %>% str_to_title()
# obj <- NormalizeData(obj, assay="RNA")

# # ## Variable features
# obj <- FindVariableFeatures(obj, assay="RNA", nfeature = 3000)

# ## cc
# obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
# obj <- CellCycleScoring(obj, g2m.features=g2m.genes, s.features=s.genes)
# obj$cc_difference <- obj$S.Score - obj$G2M.Score
# obj$Phase <- factor(obj$Phase, levels=c("G1","S","G2M"))
# table(obj$Phase, obj$Group)
# #       HSC1 HSC2 MPP1 MPP2 MPP3 MPP4 MPP5 MPP6 MPP7
# #   G1   155  176   30    3    0    5   59   78   58
# #   S      3    1   29   36   33   27   18    3   21
# #   G2M    1    0   13   33   43   46    1    0    1


# ## scale
# obj <- ScaleData(obj, assay="RNA", vars.to.regress=c('pct_counts_Mito',"cc_difference"))


# ## order
# cell_order <- c("HSC1", "MPP6", "MPP5", "MPP7", paste0("MPP", seq(4)))
# obj$Group <- factor(obj$Group, levels = cell_order)


# features <- c("Ifngr1","Ifngr2","Ifnar1","Ifnar2", "Tnfrsf1a", "Tnfrsf1b", "Fas", "Cd40", "Cd27", "Tnfrsf11a")
# DotPlot(obj, features = features, group.by = "Group", scale = TRUE) +
#     coord_flip() +
#     scale_color_gradient(low = "grey", high = "red")+
#     labs(x = "", y = "") +
#     theme_bw()+
#     theme(axis.text.x = element_text(angle = 45, hjust = 1))+
#     force_panelsizes(rows = unit(5,"cm"), cols = unit(5, "cm"))
# ggsave("plot/52.hspc.dotplot.IFN-TNFa.pdf")

# ## --------------------------------------------------------------------------------------------------





## --------------------------------------------------------------------------------------------------
## 3. Markers between groups
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










# ## boxplot with dots
# library(ggsignif)
# library(ggpubr)
# library(ggbeeswarm)

# ident.1 <- "N1"
# ident.2 <- "N2"
# sub_df <- res_ls_df[res_ls_df$Stars != "" & res_ls_df$Compare == paste0(ident.1,"-", ident.2), ] 
# selected_gene <- sub_df$gene
# data_df <- FetchData(obj %>% subset(Group %in% c(ident.1, ident.2)), vars = c(selected_gene, "Group"))
# plot_data <- tidyr::pivot_longer(data_df, cols = -Group, names_to = "Genes", values_to = "Expression")
# plot_data$Genes <- factor(plot_data$Genes, levels = selected_gene)
# comparisons <- list(c(ident.1,ident.2))

# p_ls_1 <- list()
# for (gene in selected_gene){
#     p_ls_1[[gene]] <- ggplot(plot_data[plot_data$Genes == gene,], aes(x = Group, y = Expression, color = Group)) +
#                         geom_boxplot() +
#                         scale_color_manual(values = colGroup) +
#                         geom_quasirandom(aes(color=Group),size = 3.5, alpha = 0.7)+
#                         geom_signif(comparisons = comparisons,
#                             annotations = sub_df[sub_df$gene == gene,]$Stars,  
#                             map_signif_level = TRUE, 
#                             # y_position = seq(from = 2.5, by = 0.5, length.out = length(comparisons)),  
#                             tip_length = 0.03, 
#                             size = 0.5, 
#                             textsize = 5, 
#                             color = "black") +
#                         labs(x = "", y = "Expression level") +
#                         ggtitle(gene)+
#                         theme_bw(base_size = 14)+
#                         theme(plot.title = element_text(face = "italic"),
#                               panel.grid.major = element_blank(),
#                               panel.grid.minor = element_blank(),
#                               axis.text = element_text(color = "black"),
#                               rect = element_rect(color = "black", linewidth = 1),
#                               legend.position = "none")+
#                         force_panelsizes(rows=unit(12,"cm"), cols=unit(12,"cm"))
# }

# ident.1 <- "N1"
# ident.2 <- "N3"
# sub_df <- res_ls_df[res_ls_df$Stars != "" & res_ls_df$Compare == paste0(ident.1,"-", ident.2), ] 
# selected_gene <- sub_df$gene
# data_df <- FetchData(obj %>% subset(Group %in% c(ident.1, ident.2)), vars = c(selected_gene, "Group"))
# plot_data <- tidyr::pivot_longer(data_df, cols = -Group, names_to = "Genes", values_to = "Expression")
# plot_data$Genes <- factor(plot_data$Genes, levels = selected_gene)
# comparisons <- list(c(ident.1,ident.2))

# p_ls_2 <- list()
# for (gene in selected_gene){
#     p_ls_2[[gene]] <- ggplot(plot_data[plot_data$Genes == gene,], aes(x = Group, y = Expression, color = Group)) +
#                         geom_boxplot() +
#                         scale_color_manual(values = colGroup) +
#                         geom_quasirandom(aes(color=Group),size = 3.5, alpha = 0.7)+
#                         geom_signif(comparisons = comparisons,
#                             annotations = sub_df[sub_df$gene == gene,]$Stars,  
#                             map_signif_level = TRUE, 
#                             # y_position = seq(from = 2.5, by = 0.5, length.out = length(comparisons)),  
#                             tip_length = 0.03, 
#                             size = 0.5, 
#                             textsize = 5, 
#                             color = "black") +
#                         labs(x = "", y = "Expression level") +
#                         ggtitle(gene)+
#                         theme_bw(base_size = 14)+
#                         theme(plot.title = element_text(face = "italic"),
#                               panel.grid.major = element_blank(),
#                               panel.grid.minor = element_blank(),
#                               axis.text = element_text(color = "black"),
#                               rect = element_rect(color = "black", linewidth = 1),
#                               legend.position = "none")+
#                         force_panelsizes(rows=unit(12,"cm"), cols=unit(12,"cm"))
# }

# ident.1 <- "N2"
# ident.2 <- "N3"
# sub_df <- res_ls_df[res_ls_df$Stars != "" & res_ls_df$Compare == paste0(ident.1,"-", ident.2), ] 
# selected_gene <- sub_df$gene
# data_df <- FetchData(obj %>% subset(Group %in% c(ident.1, ident.2)), vars = c(selected_gene, "Group"))
# plot_data <- tidyr::pivot_longer(data_df, cols = -Group, names_to = "Genes", values_to = "Expression")
# plot_data$Genes <- factor(plot_data$Genes, levels = selected_gene)
# comparisons <- list(c(ident.1,ident.2))

# p_ls_3 <- list()
# for (gene in selected_gene){
#     p_ls_3[[gene]] <- ggplot(plot_data[plot_data$Genes == gene,], aes(x = Group, y = Expression, color = Group)) +
#                         geom_boxplot() +
#                         scale_color_manual(values = colGroup) +
#                         geom_quasirandom(aes(color=Group),size = 3.5, alpha = 0.7)+
#                         geom_signif(comparisons = comparisons,
#                             annotations = sub_df[sub_df$gene == gene,]$Stars,  
#                             map_signif_level = TRUE, 
#                             # y_position = seq(from = 2.5, by = 0.5, length.out = length(comparisons)),  
#                             tip_length = 0.03, 
#                             size = 0.5, 
#                             textsize = 5, 
#                             color = "black") +
#                         labs(x = "", y = "Expression level") +
#                         ggtitle(gene)+
#                         theme_bw(base_size = 14)+
#                         theme(plot.title = element_text(face = "italic"),
#                               panel.grid.major = element_blank(),
#                               panel.grid.minor = element_blank(),
#                               axis.text = element_text(color = "black"),
#                               rect = element_rect(color = "black", linewidth = 1),
#                               legend.position = "none")+
#                         force_panelsizes(rows=unit(12,"cm"), cols=unit(12,"cm"))
# }
# pdf("plot/74.mpp.markers-between-grp.boxplot.pdf", height = 30, width = 50)
# print(plot_grid(plotlist = p_ls_1, ncol = 6))
# print(plot_grid(plotlist = p_ls_2, ncol = 6))
# print(plot_grid(plotlist = p_ls_3, ncol = 6))
# dev.off()




# ## heatmap
# dt <- FetchData(obj, vars = kept_gene) %>% t() %>% as.matrix()
# cell_order <- order(obj$Group)
# obj$Group[cell_order] %>% unique()
# dt <- dt[, cell_order]
# dt <- zScore(dt)

# annotation_cells <- data.frame(row.names = colnames(dt), Group = obj$Group[cell_order])
# annotation_colors <- list(Group = colGroup)

# bk <- c(seq(-2, -0.1, 0.1), seq(0, 2, 0.1))
# colHeatmap <- c(colorRampPalette(c("#5AB0C4", "#BEE2F2"))(length(bk)/2), 
#               colorRampPalette(c("#BEE2F2", "#DE482D"))(length(bk)/2))

# gaps_col <- cumsum(table(obj$Group))
# pheatmap::pheatmap(dt, cluster_row = TRUE, cluster_col = FALSE, 
#                     show_rownames = TRUE, show_colnames = FALSE, 
#                     scale = "row", 
#                     cellwidth = 400/ncol(dt), 
#                     cellheight = 300/nrow(dt), 
#                     breaks = bk, 
#                     color = colHeatmap,
#                     annotation_col = annotation_cells, 
#                     annotation_color = annotation_colors,
#                     treeheight_row = 0,
#                     gaps_col = gaps_col)


# ggplot(res_ls_df[res_ls_df$gene %in% kept_gene,], aes(x = Compare, y = gene, fill = -1*avg_log2FC)) +
#     geom_tile(na.rm = TRUE, alpha = 1) +
#     scale_fill_gradientn(colors = c("#10B3DD", "#E41A78")) +
#     geom_text(aes(label=Stars), size = 3)+
#     labs(title = paste0(""), x = "", y = "") +
#     theme(panel.grid = element_blank(), panel.background = element_blank(), axis.ticks = element_blank())+
#     # scale_x_discrete(labels = c("N2-N1", "N3-N1"))+
#     theme(axis.text = element_text(color = "black"), 
#           axis.text.x = element_text(angle = 45, hjust = 1))+
#     ggtitle("")+
#     force_panelsizes(rows=unit(8,"cm"), cols=unit(2,"cm"))
# ggsave("plot/74.mpp.markers-between-grp.tileplot.pdf", height = 10, width = 10)
## --------------------------------------------------------------------------------------------------












# ## --------------------------------------------------------------------------------------------------
# ## 3. Markers of all groups
# Idents(obj) <- "Group"

# markers <- FindAllMarkers(obj, logfc.threshold = 0, min.pct = 0)
# markers <- markers[order(-markers$avg_log2FC),]
# table(markers$cluster)

# # clean markers
# kept_genes <- rownames(markers)[!grepl("Rpl|Rps|Gm|mt-", rownames(markers))]
# markers_clean <-markers[kept_genes,]

# rds_file <- paste0("Rdata/Rds16.mpp.markers.all.rds")
# saveRDS(markers_clean, file = rds_file)
# res_ls_df <- readRDS(rds_file)
# head(res_ls_df)
# ## Stars
# res_ls_df$Stars <- cut(res_ls_df$p_val_adj,
#                            breaks = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
#                            labels = c("****", "***", "**", "*", ""))
# kept_gene <- res_ls_df[res_ls_df$Stars != "","gene"]
# res_ls_df[res_ls_df$gene == "Npm1",]


# ggplot(res_ls_df[res_ls_df$gene %in% kept_gene, ], aes(x = cluster, y = gene, fill = avg_log2FC)) +
#     geom_tile(na.rm = TRUE, alpha = 1) +
#     # scale_fill_gradientn(colors = c("#92A3D3", "#C5AFA7", "#FCB247")) +
#     scale_fill_gradient(low = "#10B3DD",  high = "#E41A78") +
#     geom_text(aes(label=Stars), size = 3)+
#     labs(title = paste0(""), x = "", y = "") +
#     theme(panel.grid = element_blank(), panel.background = element_blank(), axis.ticks = element_blank())+
#     theme(axis.text.y = element_text(color = "black"))+
#     ggtitle("")+
#     force_panelsizes(rows=unit(8,"cm"), cols=unit(1.5,"cm"))

# ## heatmap
# dt <- FetchData(obj, vars = kept_gene) %>% t() %>% as.matrix()
# cell_order <- order(obj$Group)
# obj$Group[cell_order] %>% unique()
# dt <- dt[, cell_order]
# dt <- zScore(dt)

# annotation_cells <- data.frame(row.names = colnames(dt), Group = obj$Group[cell_order])
# annotation_colors <- list(Group = colGroup)

# bk <- c(seq(-2, -0.1, 0.1), seq(0, 2, 0.1))
# colHeatmap <- c(colorRampPalette(c("#5AB0C4", "#BEE2F2"))(length(bk)/2), 
#               colorRampPalette(c("#BEE2F2", "#DE482D"))(length(bk)/2))

# gaps_col <- cumsum(table(obj$Group))
# pheatmap::pheatmap(dt, cluster_row = TRUE, cluster_col = FALSE, 
#                     show_rownames = TRUE, show_colnames = FALSE, 
#                     scale = "row", 
#                     cellwidth = 400/ncol(dt), 
#                     cellheight = 300/nrow(dt), 
#                     breaks = bk, 
#                     color = colHeatmap,
#                     annotation_col = annotation_cells, 
#                     annotation_color = annotation_colors,
#                     treeheight_row = 0,
#                     gaps_col = gaps_col)







# ## --------------------------------------------------------------------------------------------------


# feature <- "Cd27"
# # rownames(obj)[grep(feature, rownames(obj))]

# feature %in% rownames(obj)
# DimPlot(obj, group.by = "Group", reduction = "umap")/
# FeaturePlot(obj,  reduction = "umap", feature = feature)

# DotPlot(obj,  group.by = "Group", feature = feature, scale = FALSE)
