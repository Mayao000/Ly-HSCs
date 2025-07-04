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
library(ggsci) %>% suppressMessages()
library(ggalluvial) %>% suppressMessages()

### set plot options
### Define color palettes and plot themes
colLibrary <- colorRampPalette(brewer.pal(n = 11, name = "Spectral"))(11)
colGEX <- c("grey85", brewer.pal(7, "Reds"))      # color for gene expression
colCcy <- c("black", "blue", "darkorange")        # color for cellcycle phase
plotTheme <- theme_classic(base_size = 14)
# colGroup <- colorRampPalette(brewer.pal(n = 4, name = "Spectral"))(2)
colGroup <- c("#2171B5","#ED85B0", "#A1D99B")
col2rgb("#2171B5")
col2rgb("#ED85B0")

show_col(colGroup)
colTissue <- brewer.pal(4,'RdGy')[c(1,4)]
colCondition <- brewer.pal(4,'PuOr')[c(3,4)]
show_col(colCondition)
colBatch <- colorRampPalette(brewer.pal(n = 6, name = "Spectral"))(6)
show_col(colBatch)
colisPair <- c("#74C476", "#BAE4B3", "#9E9AC8")
show_col(colisPair)
# colMode <- c("#DFEBD5", "#BDB9B8")
# colMode <- c("#F96736", "#F5D206", "#0082BC")
colMode <- c("#FD3E50", "#F29837", "#7BA9C0")
show_col(colMode)
colRef <- c("black", 
            "#FF80C0", "#188386", "#808080",
            "navajowhite3", "darkorange1","forestgreen","royalblue1"
              )

## set other options after being decided
nPC <- 10
nClust <- 5
colCls <- colorRampPalette(brewer.pal(n = nClust, name = "Set1"))(nClust)
show_col(colCls)
## set python options and load modules
python_path <- '/home/mayao/.conda/envs/seurat/bin/'
use_python(python_path)
py_run_string("import numpy as np")
py_run_string("import matplotlib.pyplot as pl")


### set path
setwd("/home/mayao/LABdata/MY/pdc")

### set Seurat
options(Seurat.object.assay.version = "v5")




# ## ------------------------------------------------------------------------
# load("Rdata/02.pdc.allgrp.randomforest-predict.rdata")
# pdc$RF_predict_labels <- gsub("MPP6", "N1", pdc$RF_predict_labels)
# pdc$RF_predict_labels <- gsub("MPP5", "N2", pdc$RF_predict_labels)
# pdc$RF_predict_labels <- gsub("MPP7", "Ly-I", pdc$RF_predict_labels)
# table(pdc$Group,pdc$RF_predict_labels)

# label_order <-  c("HSC1", "N1", "N2", "Ly-I", paste0("MPP", seq(4)))
# pdc$RF_predict_labels <- factor(pdc$RF_predict_labels, levels = label_order)
# colRef <- setNames(colRef, nm = label_order)


# ## ------------------------------------------------------------------------
# ## all cells whether paired or not
# obj <- pdc
# table(obj$Group, obj$RF_predict_labels)
# pdc_res <- FetchData(obj, vars = c("Group", "RF_predict_labels"))
# pdc_res$Group <- gsub("_HSC1", "", pdc_res$Group)
# pdc_res$Group <- factor(pdc_res$Group, levels = c("Fresh", "ST", "S12"))
# table(pdc_res$Group)
# head(pdc_res)
# df_c <- ddply(pdc_res, .(Group, RF_predict_labels), summarise, Count = length(RF_predict_labels))
# df_p <- ddply(df_c, .(Group), summarise, RF_predict_labels=RF_predict_labels, Percent = Count/sum(Count)*100)
# df_p
# plot_df <- tidyr::pivot_wider(df_p, id_cols = "Group", names_from = "RF_predict_labels", values_from = "Percent")
# plot_df <- plot_df %>%
#             mutate_all(~replace(., is.na(.), 0))
# View(plot_df)
# radar_data <- rbind(rep(100, 9), rep(0, 9), plot_df)

# ## Stacked bar plot
# ggplot(df_p, aes(x = Group, y = Percent, fill = RF_predict_labels, stratum = RF_predict_labels, alluvium =  RF_predict_labels)) + 
#       geom_flow(color = "white", width = 0.5, alpha = 0.3, knot.pos=0) +
#       # geom_flow(aes(color = RF_predict_labels), fill = "white", width = 0.5, linetype = 2, alpha = 1, knot.pos=0) + 
#       # geom_flow(color = "grey", fill = "white", width = 0.5, linetype = 2, alpha = 1, knot.pos=0) + 
#       geom_col(width = 0.5, color = 'white') + 
#       scale_y_continuous(expand = c(0, 0)) +  
#       scale_fill_manual(values = colRef, name = "") +
#       # scale_color_manual(values = colRef, name = "") +
#       xlab("") + 
#       ylab("% Pairs") + 
#       theme_classic() +
#       theme(legend.title = element_blank())+
#       force_panelsizes(rows=unit(5,"cm"), cols=unit(5,"cm"))
# ggsave("plot/64.pdc.allgrp.randomforest-predict.component-bar-plot.pdf")


# # ## redar
# # library(fmsb)
# # colGroup <- setNames(colGroup, nm = c("ST", "S12", "Fresh"))
# # color_fill <- c(rgb(161, 217, 155, 70, max = 255), rgb(33, 113, 181, 70, max = 255), rgb(237, 133, 176, 70, max = 255))
# # radarchart(radar_data, 
# #            axistype = 1, 
# #            pcol = colGroup[levels(radar_data$Group)],  
# #            pfcol = color_fill ,
# #            plwd = 0.5, 
# #            cglcol = "grey", cglty = 1, 
# #            axislabcol = "black",  
# #            vlcex = 0.8 )
# # legend(x=1, y=1, legend=plot_df$Group, bty="n", pch=20, col=colGroup[levels(plot_df$Group)], text.col="black", cex=1, pt.cex=1.5)


# # col2rgb(colGroup)
# # rgb(161, 217, 155, 100, max = 255) %>% show_col

# # p_ls <- list()
# # for (grp in levels(radar_data$Group) ){
# #       p_ls[[grp]] <- radarchart(radar_data[c(1,2,which(radar_data$Group == grp)), ], axistype = 1, 
# #            pcol = colGroup[grp],  # group color
# #            pfcol = colGroup[grp],  # fill color
# #            plwd = 2,  # line width
# #            cglcol = "grey", cglty = 1,  # grid settings
# #            axislabcol = "black",  # label color
# #            vlcex = 0.8  # label size
# #       )
# # }

# # ggplot(df_p, aes(x = RF_predict_labels, y = Percent, group = Group, color = Group)) +
# #   geom_polygon(fill = NA, size = 1) +      
# #   geom_line(size = 1) +                    
# #   coord_polar() +                          
# #   theme_minimal() +
# #   labs(title = "Radar Chart using ggplot2")

# # # legend
# # legend(x = "topright", legend = rownames(radar_data)[-c(1,2)], col = c("red", "blue", "green"), lty = 1, lwd = 2)


# ## ------------------------------------------------------------------------
# ## Only pairs
# obj <- pdc %>% subset(subset=isPair=="Paired")

# # df <- FetchData(obj, vars=c("Group","Pairs","RF_predict_labels"))
# # head(df)
# # df$RF_predict_labels <- factor(df$RF_predict_labels)

# # res <- my_pair_cls_order(cls_dstr_df=df, group.by="RF_predict_labels")
# # res$Pairs_only_num <- as.character(res$Pairs_only_num)
# # res$Pair_cells <- as.character(res$Pair_cells)
# # View(res)

# # res <- res %>% arrange(Pairs, RF_predict_labels)
# # res$Pairs_only_num <- factor(res$Pairs_only_num, levels=rev(unique(res$Pairs_only_num)))
# # # s12$Facet <- rep(seq(4), each = 53)
# # # st$Facet <- rep(seq(4), times = c(46,46,46,48))

# # ggplot(res %>% subset(subset=Group=="S12"), aes(x=Pair_cells, y=Pairs_only_num, color=RF_predict_labels))+
# #             geom_point(size=2)+
# #             scale_color_manual(values = colRef[names(colRef) != "MPP3"], name = "")+
# #             # facet_grid(rows = vars(Facet))+
# #             theme_bw()+
# #             theme(axis.text.x = element_text(color = "black", size = 8), 
# #                   axis.text.y = element_text(color = "black", size = 8))+
# #             ggtitle("S12 Pairs") + 
# #             force_panelsizes(rows = unit(30, "cm"), cols = unit(3, "cm"))

# any(c(1,2) == 2)


# ## Calculate symmetric and assymetric percent
# obj$Mode <- ""
# mode_ls <- c()
# for (pair in unique(obj$Pairs)){
#       pair_rf_labels <- obj@meta.data[obj$Pairs == pair, "RF_predict_labels"]

#       if ((pair_rf_labels[1] == pair_rf_labels[2]) & all(pair_rf_labels == "HSC1")){
#             mode <- "SD"
#       }else if ((pair_rf_labels[1] != pair_rf_labels[2]) & any(pair_rf_labels == "HSC1")){
#             mode <- "AD"
#       }else if (all(pair_rf_labels != "HSC1")){
#             mode <- "BC"
#       }

#       message(pair_rf_labels, ", ", mode) 

#       mode_ls <- c(mode_ls, mode)
#       obj@meta.data[obj$Pairs == pair, "Mode"] <- mode
# }
# table(obj$Mode, obj$Group)
# obj$Mode <- factor(obj$Mode, levels = c("SD","AD","BC"))


# data <- FetchData(obj, vars = c("Group", "Mode", "RF_predict_labels"))
# data_c <- ddply(data, .(Group, Mode), summarise, Counts = length(Mode))
# data_p <- ddply(data_c, .(Group), summarise, Mode=Mode, Percent = Counts/sum(Counts)*100)

# ggplot(data_p, aes(x = "", y = Percent, fill = Mode)) +
#       geom_bar(stat = "identity", width = 1, color = "white") +
#       scale_fill_manual(values=colMode) +
#       coord_polar(theta = "y") +
#       facet_wrap(~ Group) +  
#       labs(title = "",
#             x = NULL,
#             y = NULL) +
#       geom_text(aes(label = paste0(as.character(round(Percent, digits = 1)), "%") ),             
#             position = position_stack(vjust = 0.5)) +
#       theme_void() +  # Removes background, grid, and labels
#       theme(legend.title = element_blank(), 
#             legend.text = element_text(size = 8), 
#             text = element_text(size = 12))+
#       force_panelsizes(rows = unit(5, "cm"), cols = unit(5, "cm"))
# ggsave("plot/65.pdc.allgrp.randomforest-predict.division-mode-pie.pdf")


### -----------------------------------------------------------------------------------------------------------------------
### NOT UPDATED (2503031) !!!
# ## pair predicted labels
# obj$Pair_rf_labels <- ""
# for (pair in unique(obj$Pairs)){
#       pair_rf_labels <- paste(obj@meta.data[obj$Pairs == pair, "RF_predict_labels"], collapse = "-")
#       # pair_labels_ls <- c(pair_labels_ls, pair_rf_labels)
#       obj@meta.data[obj$Pairs == pair, "Pair_rf_labels"] <- pair_rf_labels
# }
# table(obj$Pair_rf_labels)

# obj$Pair_rf_labels <- sapply(strsplit(obj$Pair_rf_labels, "-"), function(x) paste(sort(x), collapse = "-"))
# table(obj$Pair_rf_labels)
# obj$Pair_rf_labels <- factor(obj$Pair_rf_labels)

# # alluvial
# pdc_res <- FetchData(obj, vars = c("Group", "Pair_rf_labels"))[grep("HSC1", obj$Pair_rf_labels), ]
# df_c <- ddply(pdc_res, .(Group, Pair_rf_labels), summarise, Count = length(Pair_rf_labels))
# df_p <- ddply(df_c, .(Group), summarise, Pair_rf_labels=Pair_rf_labels, Percent = Count/sum(Count)*100)
# df_p
# ggplot(df_p, aes(x = Group, y = Percent, fill = Pair_rf_labels, stratum = Pair_rf_labels, alluvium =  Pair_rf_labels)) + 
#       geom_flow(width = 0.5, alpha = 0.3, knot.pos=0, color = 'white') + 
#       geom_col(width = 0.5, color = 'white') + 
#       scale_y_continuous(expand = c(0, 0)) +  
#       scale_fill_manual(values = pal_lancet("lanonc", alpha = 1)(16), name = "") +
#       xlab("") + 
#       ylab("% Pairs") + 
#       theme_classic() +
#       theme(legend.title = element_blank())+
#       force_panelsizes(rows=unit(6,"cm"), cols=unit(5,"cm"))
# ggsave("plot/66.pdc.allgrp.randomforest-predict.paired.HSC1.alluvial.pdf")


# pdc_res <- FetchData(obj, vars = c("Group", "Pair_rf_labels"))[-grep("HSC1", obj$Pair_rf_labels), ]
# df_c <- ddply(pdc_res, .(Group, Pair_rf_labels), summarise, Count = length(Pair_rf_labels))
# df_p <- ddply(df_c, .(Group), summarise, Pair_rf_labels=Pair_rf_labels, Percent = Count/sum(Count)*100)
# df_p
# ggplot(df_p, aes(x = Group, y = Percent, fill = Pair_rf_labels, stratum = Pair_rf_labels, alluvium =  Pair_rf_labels)) + 
#       geom_flow(width = 0.5, alpha = 0.3, knot.pos=0, color = 'white') + 
#       geom_col(width = 0.5, color = 'white') + 
#       scale_y_continuous(expand = c(0, 0)) +  
#       scale_fill_futurama()+
#       xlab("") + 
#       ylab("% Pairs") + 
#       theme_classic() +
#       theme(legend.title = element_blank())+
#       force_panelsizes(rows=unit(6,"cm"), cols=unit(5,"cm"))
# ggsave("plot/66.pdc.allgrp.randomforest-predict.paired.MPP.alluvial.pdf")


# ### Sankey
# df <- FetchData(obj, vars = c("Group", "Pairs", "RF_predict_labels"))
# p_ls <- list()
# for (grp in c("ST", "S12")){
#       pdc_res <- df %>% subset(subset = Group == grp)
#       pdc_pairs <- data.frame(
#                         Pair = pdc_res$Pairs[seq(1, nrow(pdc_res), 2)], # 取每对的第一个细胞
#                         Cell1 = pdc_res$RF_predict_labels[seq(1, nrow(pdc_res), 2)],
#                         Cell2 = pdc_res$RF_predict_labels[seq(2, nrow(pdc_res), 2)]
#                   )
#       p_ls[[grp]] <- ggplot(pdc_pairs, aes(axis1 = Cell1, axis2 = Cell2, y = 1)) +  # 使用 y = 1 表示所有 pair cell 的流动
#                         geom_alluvium(aes(fill = Cell1)) +
#                         geom_stratum() +
#                         scale_fill_manual(values = colRef)+
#                         geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
#                         theme_void() +
#                         labs(title = "Pair Daughter Cell Classification")+
#                         ggtitle(grp) +
#                         force_panelsizes(rows=unit(25,"cm"), cols=unit(8,"cm"))
# }
# pdf("plot/67.pdc.allgrp.randomforest-predict.paired.sankey.pdf", width = 18, height = 25)
# print(plot_grid(plotlist = p_ls, ncol = 2))
# dev.off()

# ## Details in Symmetric and Asymmetric
# pdc_res <- FetchData(obj %>% subset(subset = Mode == "Symmetric"), vars = c("Group", "Pair_rf_labels"))
# head(pdc_res )
# table(pdc_res$Group, pdc_res$Pair_rf_labels)
# df_c <- ddply(pdc_res, .(Group, Pair_rf_labels), summarise, Count = length(Pair_rf_labels))
# df_p <- ddply(df_c, .(Group), summarise, Pair_rf_labels=Pair_rf_labels, Percent = Count/sum(Count)*100)
# df_p
# ggplot(df_p, aes(x = Group, y = Percent, fill = Pair_rf_labels, stratum = Pair_rf_labels, alluvium =  Pair_rf_labels)) + 
#       geom_flow(width = 0.5, alpha = 0.3, knot.pos=0, color = 'white') + 
#       geom_col(width = 0.5, color = 'white') + 
#       scale_y_continuous(expand = c(0, 0)) +  
#       scale_fill_lancet()+
#       xlab("") + 
#       ylab("% Pairs") + 
#       theme_classic() +
#       theme(legend.title = element_blank())+
#       ggtitle("Symmetric division") +
#       force_panelsizes(rows=unit(6,"cm"), cols=unit(5,"cm"))
# ggsave("plot/68.pdc.allgrp.randomforest-predict.paired.Symmetric.alluvial.pdf")

# pdc_res <- FetchData(obj %>% subset(subset = Mode == "Asymmetric"), vars = c("Group", "Pair_rf_labels"))
# head(pdc_res )
# table(pdc_res$Group, pdc_res$Pair_rf_labels)
# df_c <- ddply(pdc_res, .(Group, Pair_rf_labels), summarise, Count = length(Pair_rf_labels))
# df_p <- ddply(df_c, .(Group), summarise, Pair_rf_labels=Pair_rf_labels, Percent = Count/sum(Count)*100)
# df_p
# ggplot(df_p, aes(x = Group, y = Percent, fill = Pair_rf_labels, stratum = Pair_rf_labels, alluvium =  Pair_rf_labels)) + 
#       geom_flow(width = 0.5, alpha = 0.3, knot.pos=0, color = 'white') + 
#       geom_col(width = 0.5, color = 'white') + 
#       scale_y_continuous(expand = c(0, 0)) +  
#       scale_fill_manual(values = c(pal_futurama("planetexpress")(12), "lightgrey"))+
#       xlab("") + 
#       ylab("% Pairs") + 
#       theme_classic() +
#       theme(legend.title = element_blank())+
#       ggtitle("Asymmetric division") +
#       force_panelsizes(rows=unit(6,"cm"), cols=unit(5,"cm"))
# ggsave("plot/68.pdc.allgrp.randomforest-predict.paired.Asymmetric.alluvial.pdf")
### -----------------------------------------------------------------------------------------------------------------------



### -----------------------------------------------------------------------------------------------------------------------
# ### Circlize Plot of AD
# library(circlize)
# all_data <- FetchData(obj, vars = c("Group", "Pairs", "Mode", "RF_predict_labels"))
# for (grp in c("ST","S12")){
#       data <- all_data %>% subset(subset = Mode == "AD" & Group == grp)

#       ## links
#       links <- data.frame(
#       from = data$RF_predict_labels[seq(1, nrow(data), by = 2)],
#       to = data$RF_predict_labels[seq(2, nrow(data), by = 2)]
#       )
#       ## sort links
#       priority <- c("HSC1", "N1", "N2", "Ly-I", "MPP1", "MPP2", "MPP3", "MPP4")
#       priority_map <- setNames(seq_along(priority), priority)
#       reorder_rows <- function(row) {
#             if (priority_map[row["from"]] > priority_map[row["to"]]) {
#                   temp <- row["from"]
#                   row["from"] <- row["to"]
#                   row["to"] <- temp
#                   }
#             return(row)
#       }
#       links_swap <- as.data.frame(t(apply(links, 1, reorder_rows)))
#       links_swap <- links_swap[order(priority_map[links_swap$from], priority_map[links_swap$to]), ]
#       unique_labels <- unique(c(links_swap$from, links_swap$to))
#       sector_order <- c("HSC1","N1", "N2", "Ly-I", "MPP1", "MPP2", "MPP3", "MPP4")

#       pdf(paste0("plot/69.pdc.",grp,".randomforest-predict.paired.AD.circos.pdf"))
#       ## initialize
#       circos.clear()
#       circos.par(gap.degree = 5)  # gap
#       circos.initialize(sectors = unique_labels, factors = unique_labels, xlim = c(0, 1))
#       ## colors
#       color_map <- c(
#       "HSC1" = "black",
#       "N1" = "#FF80C0",
#       "N2" = "#188386",
#       "Ly-I" = "#808080",
#       "MPP1" = "navajowhite3",
#       "MPP2" = "darkorange1",
#       "MPP3" = "forestgreen",
#       "MPP4" = "royalblue1"
#       )
#       links_swap$color <- color_map[links_swap$from]
#       chordDiagram(
#             links_swap[, 1:2],  
#             order = sector_order, 
#             col = links_swap$color,  
#             transparency = 0.5,
#             annotationTrack = "grid",
#             preAllocateTracks = 1
#       )
#       # 添加外圈标签并应用颜色
#       circos.trackPlotRegion(
#       ylim = c(0, 1),
#       panel.fun = function(x, y) {
#       sector.name <- get.cell.meta.data("sector.index")
#       circos.text(
#             CELL_META$xcenter, CELL_META$ycenter + 0.5, sector.name,
#             facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5),
#             col = color_map[sector.name]  
#       )
#       },
#       #   bg.col = color_map[get.cell.meta.data("sector.index")],
#       bg.border = NA  
#       )
#       dev.off()
# }



# ### Circlize Plot of BC
# library(circlize)
# all_data <- FetchData(obj, vars = c("Group", "Pairs", "Mode", "RF_predict_labels"))
# for (grp in c("S12")){
#       data <- all_data %>% subset(subset = Mode == "BC" & Group == grp)

#       ## links
#       links <- data.frame(
#       from = data$RF_predict_labels[seq(1, nrow(data), by = 2)],
#       to = data$RF_predict_labels[seq(2, nrow(data), by = 2)]
#       )
#       ## sort links
#       priority <- c("HSC1", "N1", "N2", "Ly-I", "MPP1", "MPP2", "MPP3", "MPP4")
#       priority_map <- setNames(seq_along(priority), priority)
#       reorder_rows <- function(row) {
#             if (priority_map[row["from"]] > priority_map[row["to"]]) {
#                   temp <- row["from"]
#                   row["from"] <- row["to"]
#                   row["to"] <- temp
#                   }
#             return(row)
#       }
#       links_swap <- as.data.frame(t(apply(links, 1, reorder_rows)))
#       links_swap <- links_swap[order(priority_map[links_swap$from], priority_map[links_swap$to]), ]
#       unique_labels <- unique(c(links_swap$from, links_swap$to))
#       sector_order <- c("HSC1","N1", "N2", "Ly-I", "MPP1", "MPP2", "MPP3", "MPP4")

#       pdf(paste0("plot/70.pdc.",grp,".randomforest-predict.paired.BC.circos.pdf"))
#       ## initialize
#       circos.clear()
#       circos.par(gap.degree = 5)  # gap
#       circos.initialize(sectors = unique_labels, factors = unique_labels, xlim = c(0, 1))
#       ## colors
#       color_map <- c(
#       "HSC1" = "black",
#       "N1" = "#FF80C0",
#       "N2" = "#188386",
#       "Ly-I" = "#808080",
#       "MPP1" = "navajowhite3",
#       "MPP2" = "darkorange1",
#       "MPP3" = "forestgreen",
#       "MPP4" = "royalblue1"
#       )
#       links_swap$color <- color_map[links_swap$from]
#       chordDiagram(
#             links_swap[, 1:2],  
#             order = sector_order, 
#             col = links_swap$color,  
#             transparency = 0.5,
#             annotationTrack = "grid",
#             preAllocateTracks = 1
#       )
#       # 添加外圈标签并应用颜色
#       circos.trackPlotRegion(
#       ylim = c(0, 1),
#       panel.fun = function(x, y) {
#       sector.name <- get.cell.meta.data("sector.index")
#       circos.text(
#             CELL_META$xcenter, CELL_META$ycenter + 0.5, sector.name,
#             facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5),
#             col = color_map[sector.name]  
#       )
#       },
#       #   bg.col = color_map[get.cell.meta.data("sector.index")],
#       bg.border = NA  
#       )
#       dev.off()
# }
### -----------------------------------------------------------------------------------------------------------------------






# ### ----------------------------------------------------------------------------------------------------------
# ### write out xlsx
# load("Rdata/02.pdc.allgrp.randomforest-predict.rdata")
# pdc$RF_predict_labels <- gsub("MPP6", "N1", pdc$RF_predict_labels)
# pdc$RF_predict_labels <- gsub("MPP5", "N2", pdc$RF_predict_labels)
# pdc$RF_predict_labels <- gsub("MPP7", "Ly-I", pdc$RF_predict_labels)
# table(pdc$Group,pdc$RF_predict_labels)

# label_order <-  c("HSC1", "N1", "N2", "Ly-I", paste0("MPP", seq(4)))
# pdc$RF_predict_labels <- factor(pdc$RF_predict_labels, levels = label_order)
# colRef <- setNames(colRef, nm = label_order)

# obj <- pdc
# df <- data.frame(row.names = rownames(obj@meta.data), Group = obj$Group, Pairs = obj$Pairs, isPair = obj$isPair, RF_predict_labels = obj$RF_predict_labels)
# head(df)
# df_mut <- df
# df_mut$RF_predict_pairs <- sapply(rownames(df), function(x){
#       if (df[x, "isPair"] == "Paired"){
#             pairs <- df[x, "Pairs"]
#             RF_predict_pairs <- paste( df[df$Pairs == pairs, "RF_predict_labels"], collapse = " - " )
#       }else {
#             RF_predict_pairs <- paste0(df[x, "RF_predict_labels"], " - none")
#       }
#       return(RF_predict_pairs)
# })
# summary_data_1 <- table(df_mut[df_mut$isPair == "Paired" & df_mut$Group == "ST", "RF_predict_pairs"]) %>% as.data.frame()
# summary_data_2 <- table(df_mut[df_mut$isPair == "Paired" & df_mut$Group == "S12", "RF_predict_pairs"]) %>% as.data.frame()

# xlsx_file <- paste0("12.pdc.random-forest.pair_info.xlsx")
# if ( file.exists(paste0("table/",xlsx_file)) ) {
#   file.remove( paste0("table/",xlsx_file) )
# }
# wb <- my_xlsx_create(xlsx_file)
# wb <- my_xlsx_write(wb, data = df_mut, sheet = "pdc_info", rowNames = TRUE)
# wb <- my_xlsx_write(wb, data = summary_data_1, sheet = "summary_data_ST", rowNames = FALSE)
# wb <- my_xlsx_write(wb, data = summary_data_2, sheet = "summary_data_S12", rowNames = FALSE)
# my_xlsx_save(wb, xlsx_file)

# table(obj$Group, obj$isPair)








# ## --------------------------------------------------------------------------------------------------------------------------
# library(ggsignif)
# library(ggpubr)
# library(ggbeeswarm)

# ## Check specific genes in HSC1-MPP1 & HSC1-N2 pairs
# load("Rdata/02.pdc.allgrp.randomforest-predict.rdata")
# pdc$RF_predict_labels <- gsub("MPP6", "N1", pdc$RF_predict_labels)
# pdc$RF_predict_labels <- gsub("MPP5", "N2", pdc$RF_predict_labels)
# pdc$RF_predict_labels <- gsub("MPP7", "Ly-I", pdc$RF_predict_labels)
# table(pdc$Group,pdc$RF_predict_labels)

# label_order <-  c("HSC1", "N1", "N2", "Ly-I", paste0("MPP", seq(4)))
# pdc$RF_predict_labels <- factor(pdc$RF_predict_labels, levels = label_order)
# colRef <- setNames(colRef, nm = label_order)

# obj <- pdc
# table(obj$Group, obj$RF_predict_labels)

# ## add information
# df <- data.frame(row.names = rownames(obj@meta.data), Group = obj$Group, Pairs = obj$Pairs, isPair = obj$isPair, RF_predict_labels = obj$RF_predict_labels)
# df_mut <- df
# df_mut$RF_predict_pairs <- sapply(rownames(df), function(x){
#       if (df[x, "isPair"] == "Paired"){
#             pairs <- df[x, "Pairs"]
#             RF_predict_pairs <- paste( df[df$Pairs == pairs, "RF_predict_labels"], collapse = "-" )
#       }else {
#             RF_predict_pairs <- paste0(df[x, "RF_predict_labels"], "-none")
#       }
#       return(RF_predict_pairs)
# })
# table(df_mut$RF_predict_pairs)

# ## extract HSC1-MPP1 & HSC1-N2 pairs
# sub_df <- df_mut %>% subset(RF_predict_pairs %in% c("HSC1-MPP1", "MPP1-HSC1", "HSC1-N2", "N2-HSC1"))
# sub_df$RF_predict_pairs <- gsub("MPP1-HSC1", "HSC1-MPP1", sub_df$RF_predict_pairs)
# sub_df$RF_predict_pairs <- gsub("N2-HSC1", "HSC1-N2", sub_df$RF_predict_pairs)
# table(sub_df$RF_predict_pairs)


# ## specific features
# # features <- c("Apoe", "Tgm2", "Pdzk1ip1", "Mmrn1", "Clca3a1", "Cavin2", "Atad2", "Ccnb2", "Nkg7", "Dntt", "Sell", "Tespa1")
# ## Another way to visualize difference between HSC1-MPP1 & HSC1-N2 pairs
# load("/home/mayao/LABdata/MY/hsc1_hsc2/mpp/Rdata/10.hspc-combine.rf_model_A.iml-explanation-direction.genes_res_ls.rdata")
# # features <- c(genes_res_ls[["HSC1"]]$feature,  genes_res_ls[["MPP1"]]$feature,  genes_res_ls[["N2"]]$feature) %>% unique()
# features <- Reduce(rbind, genes_res_ls)$feature %>% unique()
# features <- gsub("G_", "", features)
# features <- gsub("_", "-", features) 
# features <- features[features %in% rownames(obj)] # 125(HSC1,MPP1,N2)  All:175

# ## get data
# df <- obj[["RNA"]]$data[features, rownames(sub_df)] %>% as.matrix %>% t() %>% as.data.frame()
# df$Group <- obj@meta.data[rownames(df), "Group"]
# df$RF_predict_labels <- obj@meta.data[rownames(df), "RF_predict_labels"]

# ## reorganize data
# df_mut <- tidyr::pivot_longer(df, cols = -c(Group, RF_predict_labels), names_to = "Features", values_to = "Expression")


# ## statistical analysis - anova
# exp_data <- obj[["RNA"]]$data[features, rownames(sub_df)] %>% as.matrix
# group <- df$RF_predict_labels
# res_ls <- my_anova(expression_data = exp_data, 
#                    group = group)
# all_comparisons_df <- res_ls[["all_comparisons_df"]]
# sig_df <- res_ls[["sig_df"]]
# # sig_features <- sig_df$Gene
# sig_features <- all_comparisons_df[all_comparisons_df$Turkey_padj < 0.05, "Gene"] 
# length(unique(sig_features)) # 40

# ## Stars
# # all_comparisons_df$Stars <- cut(all_comparisons_df$BH_padj,
# #             breaks = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
# #             labels = c("****", "***", "**", "*", "ns"))
# all_comparisons_df$Stars <- cut(all_comparisons_df$Turkey_padj,
#             breaks = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
#             labels = c("****", "***", "**", "*", "ns"))

# if (length(sig_features) != 0){
#       ## boxplot
#       p_ls <- list()
#       for (gene in sig_features){
#             sub_df <- all_comparisons_df[all_comparisons_df$Gene == gene, ]
#             sub_plot_data <- df_mut[df_mut$Features == gene, ]

#             # comparisons <- strsplit(as.character(sub_df$Comparison[sub_df$BH_padj < 0.05]), "-")
#             comparisons <- strsplit(as.character(sub_df$Comparison[sub_df$Turkey_padj < 0.05]), "-")
#             comparisons <- lapply(comparisons, function(x) as.character(x))

#             if ( length(comparisons) == 0 ){
#                   comparisons <- strsplit(as.character(sub_df$Comparison), "-")
#                   comparisons <- lapply(comparisons, function(x) as.character(x))

#                   y_position <- seq(from = max(sub_plot_data$Expression)+0.5, by = 2.5, length.out = length(comparisons))
#                   # y_limit <- c(0, max(y_position)+5)
#                   y_limit <- c(0,8)

#                   p_ls[[gene]] <- ggplot(sub_plot_data, aes(x = RF_predict_labels, y = Expression, color = RF_predict_labels)) +
#                                     geom_boxplot() +
#                                     scale_color_manual(values = colRef) +
#                                     geom_signif(comparisons = comparisons,
#                                           annotations = sub_df$Stars,  
#                                           map_signif_level = TRUE, 
#                                           y_position = y_position,  
#                                           tip_length = 0.03, 
#                                           size = 0.5, 
#                                           textsize = 5, 
#                                           color = "black") +
#                                     labs(title = paste0(gene), x = "", y = "Expression level") +
#                                     scale_y_continuous(limits = y_limit)+
#                                     theme_classic()+
#                                     theme(axis.text = element_text(color = "black", size = 14), 
#                                           plot.title = element_text(color = "black", face = "italic", size = 16), 
#                                           axis.title = element_text(color = "black", size = 16),
#                                           legend.position = "none")+
#                                     force_panelsizes(rows=unit(3,"cm"), cols=unit(3.5,"cm"))

#             }else{
#                   y_position <- seq(from = max(sub_plot_data$Expression)+0.5, by = 2.5, length.out = length(comparisons))
#                   print(y_position)
#                   # y_limit <- c(0, max(y_position)+5)
#                   y_limit <- c(0,8)

#                   p_ls[[gene]] <- ggplot(sub_plot_data, aes(x = RF_predict_labels, y = Expression, color = RF_predict_labels)) +
#                                           geom_boxplot() +
#                                           scale_color_manual(values = colRef) +
#                                           geom_signif(comparisons = comparisons,
#                                                       # annotations = sub_df[sub_df$BH_padj < 0.05,]$Stars,  
#                                                       annotations = sub_df[sub_df$Turkey_padj < 0.05,]$Stars, 
#                                                       map_signif_level = TRUE, 
#                                                       y_position = y_position,  
#                                                       tip_length = 0.03, 
#                                                       size = 0.5, 
#                                                       textsize = 5, 
#                                                       color = "black") +
#                                           labs(title = paste0(gene), x = "", y = "Expression level") +
#                                           scale_y_continuous(limits = y_limit)+
#                                           theme_classic()+
#                                           theme(axis.text = element_text(color = "black", size = 14), 
#                                                       plot.title = element_text(color = "black", face = "italic", size = 16), 
#                                                       axis.title = element_text(color = "black", size = 16), 
#                                                       legend.position = "none")+
#                                           force_panelsizes(rows=unit(3,"cm"), cols=unit(3.5,"cm"))
#             }
#       }
# }

# pdf("plot/75.pdc.randomforest.iml-explained-direction-all-sig-features.HSC1-MPP1_HSC1-N2_pairs.boxplot.pdf", width = 15, height = 15)
# print(plot_grid(plotlist=p_ls, ncol = 5))
# dev.off()






# ## --------------------------------------------------------------------------------------------------------------------------
# ## Another way to visualize difference between HSC1-MPP1 & HSC1-N2 pairs
# print('start')

# my_pair_correlation <- function(seurat_obj, standard_pair, compared_pair, features, group_by){
#       sub_obj <- seurat_obj %>% subset(Correlation_pair_type %in% c(standard_pair, compared_pair))
#       print(table(sub_obj$Correlation_pair_type))

#       ## statistical analysis using Seurat
#       Idents(sub_obj) <- group_by
#       deg_ls <- list()
#       for ( temp in c(standard_pair, compared_pair) ){

#             temp_obj <- sub_obj %>% subset(Correlation_pair_type == temp) 
#             ident.1 <- str_split(temp, "[-]", simplify = TRUE)[,1]
#             ident.2 <- str_split(temp, "[-]", simplify = TRUE)[,2]
#             if (ident.1 == ident.2){
#                   ident.1 <- paste0(ident.1, ".1")
#                   ident.2 <- paste0(ident.1, ".2")
#             }
#             message("ident.1 is ", ident.1, "\nident.2 is ", ident.2, "\n")

#             deg_ls[[temp]] <- FindMarkers(temp_obj, 
#                                     recorrect_umi=FALSE, 
#                                     ident.1 = ident.1,
#                                     ident.2 = ident.2,
#                                     features = features,
#                                     logfc.threshold = 0, 
#                                     min.pct = 0, 
#                                     return.thresh = 1,
#                                     min.cells.feature = 0)
#             deg_ls[[temp]]$Cluster <- ifelse(deg_ls[[temp]]$avg_log2FC > 0, ident.1, ident.2)
#             deg_ls[[temp]] <- deg_ls[[temp]][ order( deg_ls[[temp]]$avg_log2FC, decreasing = TRUE), ]
#       }


#       gene_order <- deg_ls[[standard_pair]] %>% rownames()

#       df_long <- data.frame(
#       GeneIndex = 1:length(gene_order),
#             standard = deg_ls[[standard_pair]][gene_order, ]$avg_log2FC,
#             compared = deg_ls[[compared_pair]][gene_order, ]$avg_log2FC
#       ) %>% tidyr::pivot_longer(cols = -GeneIndex, names_to = "Group", values_to = "log2FC")

#       p <- ggplot(df_long, aes(x = GeneIndex, y = log2FC, color = Group, fill = Group)) +
#                   geom_point(size = 0.5)+
#                   geom_smooth(se = TRUE, span = 0.2, alpha = 0.1) +  
#                   scale_color_manual(values = c("#0094D6", "#F7AE24"))+
#                   scale_fill_manual(values =c("#0094D6", "#F7AE24"))+
#                   labs(x = "Gene index \n (ordered by HSC1-MPP1 log2FC)", y = "Absolute log2(Fold change)",
#                         title = "") +
#                   theme_classic()+
#                   theme(axis.text = element_text(size = 14, color = "black"), 
#                         axis.title = element_text(size = 14, color = "black"),
#                         plot.title = element_text(size = 14, color = "black"))+
#                   force_panelsizes(rows = unit(5, 'cm'), cols = unit(7, "cm"))

#       vec1 <- deg_ls[[standard_pair]][gene_order, ]$avg_log2FC
#       vec2 <- deg_ls[[compared_pair]][gene_order, ]$avg_log2FC
#       cor_obs <- cor(vec1, vec2, method = "spearman") # 0.039

#       res_ls <- c(list(p), list(cor_obs))
#       names(res_ls) <- c("p", "cor_obs")
#       return(res_ls)
# }


# ## Another way to visualize difference between HSC1-MPP1 & HSC1-N2 pairs
# load("/home/mayao/LABdata/MY/hsc1_hsc2/mpp/Rdata/10.hspc-combine.rf_model_A.iml-explanation-direction.genes_res_ls.rdata")
# # features <- c(genes_res_ls[["HSC1"]]$feature,  genes_res_ls[["MPP1"]]$feature,  genes_res_ls[["N2"]]$feature) %>% unique()
# features <- Reduce(rbind, genes_res_ls)$feature %>% unique()
# features <- gsub("G_", "", features)
# features <- gsub("_", "-", features) 
# features <- features[features %in% rownames(obj)] # 125(HSC1,MPP1,N2)  All:175


# load("Rdata/02.pdc.allgrp.randomforest-predict.rdata")
# pdc$RF_predict_labels <- gsub("MPP6", "N1", pdc$RF_predict_labels)
# pdc$RF_predict_labels <- gsub("MPP5", "N2", pdc$RF_predict_labels)
# pdc$RF_predict_labels <- gsub("MPP7", "LyI", pdc$RF_predict_labels)
# table(pdc$Group,pdc$RF_predict_labels)

# label_order <-  c("HSC1", "N1", "N2", "LyI", paste0("MPP", seq(4)))
# pdc$RF_predict_labels <- factor(pdc$RF_predict_labels, levels = label_order)
# colRef <- setNames(colRef, nm = label_order)

# obj <- pdc
# table(obj$Group, obj$RF_predict_labels)

# ## add information
# df <- data.frame(row.names = rownames(obj@meta.data), Group = obj$Group, Pairs = obj$Pairs, isPair = obj$isPair, RF_predict_labels = obj$RF_predict_labels)
# obj$RF_predict_pairs <- sapply(rownames(df), function(x){
#       if (df[x, "isPair"] == "Paired"){
#             pairs <- df[x, "Pairs"]
#             RF_predict_pairs <- paste( df[df$Pairs == pairs, "RF_predict_labels"], collapse = "-" )
#       }else {
#             RF_predict_pairs <- paste0(df[x, "RF_predict_labels"], "-none")
#       }
#       return(RF_predict_pairs)
# })
# table(obj$RF_predict_pairs)


# paired_obj <- obj %>% subset(isPair == "Paired")
# paired_obj$RF_predict_labels_new <- as.character(paired_obj$RF_predict_labels)

# meta <- paired_obj@meta.data %>% as.data.frame()
# paired_obj$RF_predict_labels_new[seq(2, nrow(meta), 2)] <- sapply(seq(2, nrow(meta), 2), function(x) {
#       temp <- meta[x, "RF_predict_pairs"]
#       ident.1 <- str_split(temp, "[-]", simplify = TRUE)[,1]
#       ident.2 <- str_split(temp, "[-]", simplify = TRUE)[,2]

#       if (ident.1 == ident.2){
#             RF_predict_labels_new <- paste0(ident.2, ".2")
#       }else{
#             RF_predict_labels_new <- ident.2
#       }
      
#       return(RF_predict_labels_new)
# })
# table(paired_obj$RF_predict_labels_new)

# pairs <- paste0( paired_obj$RF_predict_labels_new[ seq(1, nrow(meta), 2) ], "-", paired_obj$RF_predict_labels_new[ seq(2, nrow(meta), 2) ] )
# paired_obj$RF_predict_pairs_new <- rep(pairs, each = 2)
# paired_obj$Correlation_pair_type <- sapply(paired_obj$RF_predict_pairs_new, function(x) {
#   sorted <- sort(unlist(strsplit(x, "-")))
#   paste(sorted, collapse = "-")
# })
# table(paired_obj$Correlation_pair_type)




# group_by <- "RF_predict_labels_new"
# pair_pool <- names( which(table(paired_obj$Correlation_pair_type) > 6 ) )
# cor_obs_ls <- list()
# standard_pair_ls <- c()
# for (standard_pair in pair_pool){
#       standard_pair_ls <- c(standard_pair_ls, standard_pair)
#       compared_pair_ls <- pair_pool[ -which(pair_pool %in% standard_pair_ls) ]

#       for(compared_pair in compared_pair_ls){
#             res_ls <- my_pair_correlation(seurat_obj = paired_obj, standard_pair = standard_pair, compared_pair = compared_pair, features = features, group_by = group_by)
#             cor_obs_ls[[paste0(standard_pair,"_", compared_pair)]] <- res_ls[["cor_obs"]]
#       }
      
# }
# cor_obs_ls


# ## get null distribution from HSC1-HSC1 pairs using permutation test
# my_permutation <- function(seurat_obj, features, permutation_num = 1000, n_cell){
#       cor_res_ls <- c()
#       for (num in 1:permutation_num){
#             set.seed(num)
#             message("Test: ", num)
#             labels_1 <- sample(rep(c("PDC1", "PDC2"), each = n_cell/2))
#             # labels_1 <- character(n_cell)
#             # for (i in seq(1, n, by = 2)) {
#             #       pair_label <- sample(c("PDC1", "PDC2"))
#             #       labels_1[i] <- pair_label[1]      
#             #       labels_1[i + 1] <- pair_label[2]  
#             # }

#             seurat_obj$PDC_label <- as.character(seurat_obj$RF_predict_labels)
#             seurat_obj@meta.data[seurat_obj$isPair == "Paired", ]$PDC_label <- labels_1

#             set.seed(num+1234)
#             labels_2 <- sample(rep(c("PDC1", "PDC2"), each = n_cell/2))
#             # labels_2 <- character(n_cell)
#             # for (i in seq(1, n, by = 2)) {
#             #       pair_label <- sample(c("PDC1", "PDC2"))
#             #       labels_2[i] <- pair_label[1]      
#             #       labels_2[i + 1] <- pair_label[2]  
#             # }

#             seurat_obj$PDC_label_2 <- as.character(seurat_obj$RF_predict_labels)
#             seurat_obj@meta.data[seurat_obj$isPair == "Paired", ]$PDC_label_2 <- labels_2

#             ident.1 <- "PDC1"
#             ident.2 <- "PDC2"
#             deg_ls <- list()
#             for (group_by in c("PDC_label", "PDC_label_2")){
#                   Idents(seurat_obj) <- group_by
#                   temp_obj <- seurat_obj %>% subset(isPair == "Paired") 
#                   deg_ls[[group_by]] <- FindMarkers(temp_obj, 
#                                     recorrect_umi=FALSE, 
#                                     ident.1 = ident.1,
#                                     ident.2 = ident.2,
#                                     features = features,
#                                     logfc.threshold = 0, 
#                                     min.pct = 0, 
#                                     return.thresh = 1,
#                                     min.cells.feature = 0, 
#                                     verbose = FALSE)
#                   deg_ls[[group_by]]$Cluster <- ifelse(deg_ls[[group_by]]$avg_log2FC > 0, ident.1, ident.2)
#                   # deg_ls[[group_by]] <- deg_ls[[group_by]][ order( abs(deg_ls[[group_by]]$avg_log2FC), decreasing = TRUE), ]
#                   deg_ls[[group_by]] <- deg_ls[[group_by]][ order( deg_ls[[group_by]]$avg_log2FC, decreasing = TRUE), ]
#             }     

#             gene_order <- rownames( deg_ls[["PDC_label"]] )
            
#             # vec1 <- abs(deg_ls[["PDC_label"]][gene_order, ]$avg_log2FC)
#             # vec2 <- abs(deg_ls[["PDC_label_2"]][gene_order, ]$avg_log2FC)

#             vec1 <- deg_ls[["PDC_label"]][gene_order, ]$avg_log2FC
#             vec2 <- deg_ls[["PDC_label_2"]][gene_order, ]$avg_log2FC

#             cor_res <- cor( vec1, vec2, method = "spearman")
#             cor_res_ls <- c(cor_res_ls, cor_res)
#       }
      
#       print(cor_res_ls)
#       null_val <- quantile(cor_res_ls, probs = c(0.025, 0.975))
#       res_ls <- c(list(cor_res_ls), list(null_val))
#       names(res_ls) <- c("cor_res_ls", "null_val")
#       return(res_ls)
# }

# res_ls <- my_permutation(seurat_obj=obj, features = features, permutation_num = 1000, n_cell = 398)
# # saveRDS(res_ls, file = "Rdata/04.pdc.randomforest-pos-correlation-distribution.rds")
# # saveRDS(res_ls, file = "Rdata/05.pdc.randomforest-null-distribution.HSC1MPP1N2-top100-features.rds")
# saveRDS(res_ls, file = "Rdata/06.pdc.randomforest-null-distribution.all-top100-features.rds")



# res_ls <- readRDS("Rdata/06.pdc.randomforest-null-distribution.all-top100-features.rds")
# cor_res_ls <- res_ls[["cor_res_ls"]]
# null_val <- res_ls[["null_val"]]
# null_val
# #  2.5%      97.5% 
# # -0.2649790  0.2664447

# ## all top100 features
# #  2.5%      97.5% 
# # -0.2315236  0.2548558


# cor_obs_pls <- lapply(cor_obs_ls, function(x){
#       pval_two_tail <- mean(abs(cor_res_ls - mean(cor_res_ls)) >= abs(x - mean(cor_res_ls)))
# }) 
# cor_obs_pls[ cor_obs_pls < 0.05 ] # 0.009
# cor_obs_ls["HSC1-MPP1_HSC1-N2"] # 0.3373394

# # $`HSC1-MPP1_HSC1-N2`
# # [1] 0.014


# pair_ls <- my_pair_correlation(seurat_obj = paired_obj, standard_pair = "HSC1-MPP1", compared_pair = "HSC1-N2", features = features, group_by = group_by)
# p1 <- pair_ls[["p"]] +
#             scale_color_manual(values = c("#F7AE24","#0094D6"), labels = c("HSC1-N2","HSC1-MPP1"))+
#             scale_fill_manual(values = c("#F7AE24","#0094D6"), labels = c("HSC1-N2","HSC1-MPP1"))

# cor_obs <- cor_obs_ls[["HSC1-MPP1_HSC1-N2"]] 
# p2 <- ggplot(data.frame(x = cor_res_ls), aes(x = x)) +
#             geom_histogram(aes(y = ..density..), bins = 30, fill = "lightblue", color = "black", alpha = 0.6) +
#             geom_density(color = "blue", size = 0.5) +
#             geom_vline(xintercept = cor_obs, color = "red", size = 0.5, linetype = 2) +
#             geom_vline(xintercept = null_val, color = "blue", size = 0.5, linetype = 2) +
#             scale_x_continuous(limits = c(-1, 1))+
#             labs(title = "Null distribution", x = "Spearman correlation", y = "Density") +
#             theme_classic()+
#             theme(axis.text = element_text(size = 14, color = "black"), 
#                   axis.title = element_text(size = 14, color = "black"),
#                   plot.title = element_text(size = 14, color = "black"))+
#             force_panelsizes(rows = unit(5, 'cm'), cols = unit(7, "cm"))

# pdf("plot/77.pdc.randomforest.iml-explained-all-features.HSC1-MPP1_HSC1-N2.log2fc-distribution.null-distribution.pdf", width = 15, height = 10)
# print(plot_grid(plotlist=c(list(p1), list(p2)), ncol = 2))
# dev.off()




# ### 
# ## get null distribution from HSC1-HSC1 pairs using permutation test
# seurat_obj
# my_permutation_2 <- function(seurat_obj, features, grp_1 = "HSC1", grp_2 = "MPP1", grp_3 = "N2", permutation_num = 1000){
#       cell_names_1 <- seurat_obj@meta.data[seurat_obj$RF_predict_labels == grp_1, ] %>% rownames()
#       cell_names_2 <- seurat_obj@meta.data[seurat_obj$RF_predict_labels == grp_2, ] %>% rownames()
#       cell_names_3 <- seurat_obj@meta.data[seurat_obj$RF_predict_labels == grp_3, ] %>% rownames()

#       seurat_obj@meta.data[grp_2_cells, "RF_predict_labels"]
#       c(seurat_obj@meta.data[grp_1_cells, "RF_predict_labels"], seurat_obj@meta.data[grp_2_cells, "RF_predict_labels"])
#       temp_obj$RF_predict_labels

#       cor_res_ls <- c()
#       for (num in 1:permutation_num){
#             set.seed(num)
#             message("Test: ", num)
#             grp_1_cells <-  sample(cell_names_1, size = 15)
#             grp_2_cells <-  sample(cell_names_2, size = 15)

#             temp_obj <- seurat_obj[, c(grp_1_cells, grp_2_cells)]
#             temp_obj$PDC_label <- ifelse(temp_obj$RF_predict_labels == grp_1, "PDC1", "PDC2")
#             Idents(temp_obj) <- "PDC_label"


#             set.seed(num+1234)
#             grp_1_cells <-  sample(cell_names_1, size = 14)
#             grp_3_cells <-  sample(cell_names_3, size = 14)

#             temp_obj_2 <- seurat_obj[, c(grp_1_cells, grp_3_cells)]
#             temp_obj_2$PDC_label <- ifelse(temp_obj_2$RF_predict_labels == grp_1, "PDC1", "PDC2")
#             Idents(temp_obj_2) <- "PDC_label"
            

#             ident.1 <- "PDC1"
#             ident.2 <- "PDC2"
#             deg_ls[[paste0(grp_1, "-", grp_2)]] <- FindMarkers(temp_obj, 
#                                                       recorrect_umi=FALSE, 
#                                                       ident.1 = ident.1,
#                                                       ident.2 = ident.2,
#                                                       features = features,
#                                                       logfc.threshold = 0, 
#                                                       min.pct = 0, 
#                                                       return.thresh = 1,
#                                                       min.cells.feature = 0, 
#                                                       verbose = FALSE)
#             deg_ls[[paste0(grp_1, "-", grp_3)]] <- FindMarkers(temp_obj_2, 
#                                                       recorrect_umi=FALSE, 
#                                                       ident.1 = ident.1,
#                                                       ident.2 = ident.2,
#                                                       features = features,
#                                                       logfc.threshold = 0, 
#                                                       min.pct = 0, 
#                                                       return.thresh = 1,
#                                                       min.cells.feature = 0, 
#                                                       verbose = FALSE)
#             deg_ls[[paste0(grp_1, "-", grp_2)]] <- deg_ls[[paste0(grp_1, "-", grp_2)]][ order( deg_ls[[paste0(grp_1, "-", grp_2)]]$avg_log2FC, decreasing = TRUE), ]
#             gene_order <- rownames( deg_ls[[paste0(grp_1, "-", grp_2)]] )

#             vec1 <- deg_ls[[paste0(grp_1, "-", grp_2)]][gene_order, ]$avg_log2FC
#             vec2 <- deg_ls[[paste0(grp_1, "-", grp_3)]][gene_order, ]$avg_log2FC

#             cor_res <- cor( vec1, vec2, method = "spearman")
#             cor_res_ls <- c(cor_res_ls, cor_res)
#       }
#       print(cor_res_ls)
#       null_val <- quantile(cor_res_ls, probs = c(0.025, 0.975))
#       res_ls <- c(list(cor_res_ls), list(null_val))
#       names(res_ls) <- c("cor_res_ls", "null_val")
#       return(res_ls)
# }

# res_ls_2 <- my_permutation_2(seurat_obj = obj, features = features, grp_1 = "HSC1", grp_2 = "MPP1", grp_3 = "N2", permutation_num = 1000)
# # res_ls <- readRDS("Rdata/06.pdc.randomforest-null-distribution.all-top100-features.rds")
# cor_res_ls <- res_ls_2[["cor_res_ls"]]
# null_val <- res_ls_2[["null_val"]]
# null_val


# pval_two_tail <- mean(abs(cor_res_ls - mean(cor_res_ls)) >= abs(cor_obs - mean(cor_res_ls)))






## -----------------------------------------------------------------------------------------------
## ss of log2FC between each pair of cells
print("start")
# my_pairwise_ss_nulldistribution <- function(seu, gene_list, pair_num = 15, permutation_num = 1000) {
#       # seu: Seurat object
#       # gene_list: vector of gene names (length ~20)
#       # cell_pairs: a 15x2 matrix or data.frame of cell barcodes (each row is a pair)
      
#       ss_vec_ls <- list()
#       ss_sum_ls <- c()
#       for (test_n in 1:permutation_num){
#             set.seed(test_n)
#             message("Test ", test_n)
      
#             # pairs
#             pair_name <- paste0("Pair_", seq(pair_num))

#             # get meta
#             meta <- seu@meta.data %>% as.data.frame()

#             # get log-normalized expression matrix
#             expr_mat <- GetAssayData(seu, assay = "RNA", slot = "data")  # log1p(normalized)

#             # initialize matrix
#             log2FC_matrix <- matrix(NA, nrow = length(gene_list), ncol = pair_num)
#             rownames(log2FC_matrix) <- gene_list
#             colnames(log2FC_matrix) <- pair_name

#             # randomly sample cells
#             cells <- sample(colnames(seu), size = pair_num*2)
#             cell_pair_df <- data.frame(row.names = cells, pair_name = rep(pair_name, each = 2), pair_ident = rep(c("PDC1", "PDC2"), pair_num))


#             # calculate log2FC
#             for (g in seq_along(gene_list)) {
#                   gene <- gene_list[g]
#                   for (i in pair_name) {
#                         c1 <- cell_pair_df[cell_pair_df$pair_name == i & cell_pair_df$pair_ident == "PDC1", ] %>% rownames()
#                         c2 <- cell_pair_df[cell_pair_df$pair_name == i & cell_pair_df$pair_ident == "PDC2", ] %>% rownames()

#                         if (gene %in% rownames(expr_mat) && c1 %in% colnames(expr_mat) && c2 %in% colnames(expr_mat)) {
#                               log2FC_matrix[g, i] <- expr_mat[gene, c1] - expr_mat[gene, c2]  
#                         } else {
#                               warning(paste("Missing gene or cell:", gene, c1, c2))
#                         }
                  
#                   }
#             }
#             ss_vector <- rowSums(log2FC_matrix, na.rm = TRUE) 
#             ss_vec_ls[[paste0("Test_", test_n)]] <- ss_vector

#             ss_sum <- ss_vector %>% abs() %>% sum()
#             ss_sum_ls <- c(ss_sum_ls, sum(ss_sum))
#       }
#       ss_vec_df <- Reduce(rbind, ss_vec_ls) %>% t() %>% as.data.frame()
#       rownames(ss_vec_df) <- gene_list
#       colnames(ss_vec_df) <- paste0("null_curve_", seq(permutation_num))
#       print(head(ss_vec_df))
#       return(list(null_ss_df = ss_vec_df, null_ss_sum = ss_sum_ls))
# }


my_pairwise_ss_nulldistribution <- function(seu, genes_to_sample, gene_list_size, cell_pairs, ident.1, ident.2, permutation_num = 1000) {
      # seu: Seurat object
      # gene_list: vector of gene names (length ~20)
      # cell_pairs: a 15x2 matrix or data.frame of cell barcodes (each row is a pair)

      ss_vec_ls <- list()
      ss_sum_ls <- c()
      for (test_n in 1:permutation_num){
            # seu: Seurat object
            # gene_list: vector of gene names (length ~20)
            # cell_pairs: a 15x2 matrix or data.frame of cell barcodes (each row is a pair)

            set.seed(test_n)
            message("Test ", test_n)
            # get meta
            meta <- seu@meta.data %>% as.data.frame()

            # get log-normalized expression matrix
            expr_mat <- GetAssayData(seu, assay = "RNA", slot = "data")  # log1p(normalized)

            # sample genes
            gene_list <- sample(genes_to_sample, size = gene_list_size)

            # initialize matrix
            log2FC_matrix <- matrix(NA, nrow = length(gene_list), ncol = length(cell_pairs))
            rownames(log2FC_matrix) <- gene_list
            colnames(log2FC_matrix) <- cell_pairs

            # calculate log2FC
            for (g in seq_along(gene_list)) {
                  gene <- gene_list[g]
                  for (i in cell_pairs) {
                        c1 <- meta[meta$Pairs == i & meta$RF_predict_labels == ident.1, ] %>% rownames()
                        c2 <- meta[meta$Pairs == i & meta$RF_predict_labels == ident.2, ] %>% rownames()

                        if (gene %in% rownames(expr_mat) && c1 %in% colnames(expr_mat) && c2 %in% colnames(expr_mat)) {
                              log2FC_matrix[g, i] <- expr_mat[gene, c1] - expr_mat[gene, c2]  
                        } else {
                              warning(paste("Missing gene or cell:", gene, c1, c2))
                        }
                  
                  }
            }
            # ss_vector <- rowSums(log2FC_matrix, na.rm = TRUE)
            ss_vector <- rowMeans(log2FC_matrix, na.rm = TRUE)
            ss_vector <- ss_vector[order(ss_vector, decreasing =TRUE)]
            ss_vec_ls[[paste0("Test_", test_n)]] <- ss_vector

            ss_sum <- ss_vector %>% abs() %>% sum()
            ss_sum_ls <- c(ss_sum_ls, sum(ss_sum))
      }
      ss_vec_df <- Reduce(rbind, ss_vec_ls) %>% t() %>% as.data.frame()
      rownames(ss_vec_df) <- seq(nrow(ss_vec_df))
      colnames(ss_vec_df) <- paste0("null_curve_", seq(permutation_num))
      print(head(ss_vec_df))
      return(list(null_ss_df = ss_vec_df, null_ss_sum = ss_sum_ls))
}



## ss of log2FC between each pair of cells
my_compute_pairwise_log2FC_ss <- function(seu, gene_list, cell_pairs, ident.1, ident.2) {
      # seu: Seurat object
      # gene_list: vector of gene names (length ~20)
      # cell_pairs: a 15x2 matrix or data.frame of cell barcodes (each row is a pair)

      # get meta
      meta <- seu@meta.data %>% as.data.frame()

      # get log-normalized expression matrix
      expr_mat <- GetAssayData(seu, assay = "RNA", slot = "data")  # log1p(normalized)

      # initialize matrix
      log2FC_matrix <- matrix(NA, nrow = length(gene_list), ncol = length(cell_pairs))
      rownames(log2FC_matrix) <- gene_list
      colnames(log2FC_matrix) <- cell_pairs

      # calculate log2FC
      for (g in seq_along(gene_list)) {
            gene <- gene_list[g]
            for (i in cell_pairs) {
                  c1 <- meta[meta$Pairs == i & meta$RF_predict_labels == ident.1, ] %>% rownames()
                  c2 <- meta[meta$Pairs == i & meta$RF_predict_labels == ident.2, ] %>% rownames()
                  # message(meta[c1, "Pairs"], ", ", meta[c1, "RF_predict_labels"], "\n")
                  # message(meta[c2, "Pairs"], ", ", meta[c2, "RF_predict_labels"], "\n")

                  if (gene %in% rownames(expr_mat) && c1 %in% colnames(expr_mat) && c2 %in% colnames(expr_mat)) {
                        log2FC_matrix[g, i] <- expr_mat[gene, c1] - expr_mat[gene, c2]  
                  } else {
                        warning(paste("Missing gene or cell:", gene, c1, c2))
                  }
            
            }
      }
      # ss_vector <- rowSums(log2FC_matrix, na.rm = TRUE)
      ss_vector <- rowMeans(log2FC_matrix, na.rm = TRUE)
      ss_df <- data.frame(row.names = rownames(log2FC_matrix), features = rownames(log2FC_matrix), ss_vector = ss_vector)
      ss_df <- ss_df[order(ss_df$ss_vector, decreasing =TRUE), ]

      ss_sum <- ss_vector  %>% abs() %>% sum()
      return(list(ss_df = ss_df, ss_sum = ss_sum))
}


compare_curve_no_align <- function(real_curve, null_curves, method = c("tss"), alternative = "two.sided") {
      method <- match.arg(method)
      n_null <- ncol(null_curves)

      # # define method function
      # cosine_sim <- function(x, y) sum(x * y) / (sqrt(sum(x^2)) * sqrt(sum(y^2)))
      # cor_sim <- function(x, y) cor(x, y, method = "spearman")
      tss_cal <- function(x,y) sum((x-y)^2)

      # define ref_curve
      # ref_curve <- rowMeans(null_curves)  
      ref_curve <- rep(0, length(real_curve))

      # distance between null_curve and ref_curve
      null_scores <- numeric(n_null)
      for (i in 1:n_null) {
            null_curve <- null_curves[, i]
            null_scores[i] <- switch(method,tss = tss_cal(null_curve, ref_curve))
                                    # cosine = cosine_sim(null_curve, ref_curve),
                                    # correlation = cor_sim(null_curve, ref_curve))
      }

      # distance between real_curve and ref_curve
      obs_score <- switch(method,tss = tss_cal(real_curve, ref_curve))

      # p value
      # Cosine / cor: smaller, more similar → left tail/two sided
      if (alternative == "two.sided") {
            pval <- mean(abs(null_scores) >= abs(obs_score))
      } else if (alternative == "less") {
            pval <- mean(null_scores <= obs_score)
      } else {
            pval <- mean(null_scores >= obs_score)
      }

      return(list(method = method, obs_score = obs_score, null_scores = null_scores, p_value = pval))
}


load("Rdata/02.pdc.allgrp.randomforest-predict.rdata")
pdc$RF_predict_labels <- gsub("MPP6", "N1", pdc$RF_predict_labels)
pdc$RF_predict_labels <- gsub("MPP5", "N2", pdc$RF_predict_labels)
pdc$RF_predict_labels <- gsub("MPP7", "LyI", pdc$RF_predict_labels)
table(pdc$Group,pdc$RF_predict_labels)

label_order <-  c("HSC1", "N1", "N2", "LyI", paste0("MPP", seq(4)))
pdc$RF_predict_labels <- factor(pdc$RF_predict_labels, levels = label_order)
colRef <- setNames(colRef, nm = label_order)

obj <- pdc
table(obj$Group, obj$RF_predict_labels)

## add information
df <- data.frame(row.names = rownames(obj@meta.data), Group = obj$Group, Pairs = obj$Pairs, isPair = obj$isPair, RF_predict_labels = obj$RF_predict_labels)
obj$RF_predict_pairs <- sapply(rownames(df), function(x){
      if (df[x, "isPair"] == "Paired"){
            pairs <- df[x, "Pairs"]
            RF_predict_pairs <- paste( df[df$Pairs == pairs, "RF_predict_labels"], collapse = "-" )
      }else {
            RF_predict_pairs <- paste0(df[x, "RF_predict_labels"], "-none")
      }
      return(RF_predict_pairs)
})
table(obj$RF_predict_pairs)


## change pair order to make them unique
RF_predict_pairs <- sapply(obj$RF_predict_pairs, function(x) {
  sorted <- sort(unlist(strsplit(x, "-")))
  paste(sorted, collapse = "-")
})
table(RF_predict_pairs)

obj$RF_predict_pairs <- RF_predict_pairs


## Features
load("/home/mayao/LABdata/MY/hsc1_hsc2/mpp/Rdata/10.hspc-combine.rf_model_A.iml-explanation-direction.genes_res_ls.rdata")
features <- c(genes_res_ls[["HSC1"]]$feature,  genes_res_ls[["MPP1"]]$feature,  genes_res_ls[["N2"]]$feature) %>% unique()
# features <- Reduce(rbind, genes_res_ls)$feature %>% unique()
features <- gsub("G_", "", features)
features <- gsub("_", "-", features) 
features <- features[features %in% rownames(obj)] # 125(HSC1,MPP1,N2)  All:175

## genes_to_sample
genes_to_sample <- readRDS("/home/mayao/LABdata/MY/hsc1_hsc2/mpp/Rdata/04.hspc-combine.features.randomforest.rds")
genes_to_sample <- genes_to_sample[genes_to_sample %in% rownames(obj)] # 1321

## get ss_df
## get HSC1-MPP1 & HSC1-N2 pairs
ss_df_ls <- list()
null_ss_df_ls <- list()
# for (pair_type in c("HSC1-MPP1", "HSC1-N2")){
#       message(pair_type)
#       paired_obj <- obj %>% subset(RF_predict_pairs == pair_type)
#       ident.1 <- str_split(pair_type, "[-]", simplify = TRUE)[,1]
#       ident.2 <- str_split(pair_type, "[-]", simplify = TRUE)[,2]
#       ss_res_ls <- my_compute_pairwise_log2FC_ss(seu = paired_obj, gene_list = features, 
#                                                 cell_pairs = unique(paired_obj$Pairs),   
#                                                 ident.1 = ident.1, ident.2 = ident.2)
#       ss_df <- ss_res_ls[["ss_df"]]
#       ss_df_ls[[pair_type]] <- ss_df

#       gene_order <- rownames(ss_df)
#       null_res <- my_pairwise_ss_nulldistribution(seu=obj, gene_list=gene_order, pair_num = 15, permutation_num = 1000)
#       saveRDS(null_res, file = paste0("Rdata/10.pdc.randomforest.ss-null-distribution.125-hsc1mpp1n2-features.",pair_type,"_ordered.rds"))

#       null_ss_df <- null_res[["null_ss_df"]]
#       null_ss_df_ls[[pair_type]] <- null_ss_df
# }


for (pair_type in c("HSC1-MPP1", "HSC1-N2")){
      message(pair_type)
      paired_obj <- obj %>% subset(RF_predict_pairs == pair_type)
      ident.1 <- str_split(pair_type, "[-]", simplify = TRUE)[,1]
      ident.2 <- str_split(pair_type, "[-]", simplify = TRUE)[,2]
      ss_res_ls <- my_compute_pairwise_log2FC_ss(seu = paired_obj, gene_list = features, 
                                                cell_pairs = unique(paired_obj$Pairs),   
                                                ident.1 = ident.1, ident.2 = ident.2)
      ss_df <- ss_res_ls[["ss_df"]]
      ss_df_ls[[pair_type]] <- ss_df

      gene_order <- rownames(ss_df)
      null_res <- my_pairwise_ss_nulldistribution(seu = paired_obj, genes_to_sample = genes_to_sample, 
                                                  gene_list_size = nrow(ss_df), 
                                                  cell_pairs = unique(paired_obj$Pairs), 
                                                  ident.1 = ident.1, ident.2 = ident.2, 
                                                  permutation_num = 1000)
      saveRDS(null_res, file = paste0("Rdata/12.pdc.randomforest.rowMeans-ss.ss-null-distribution.random-features.",pair_type,"_ordered.rds"))

      null_ss_df <- null_res[["null_ss_df"]]
      null_ss_df_ls[[pair_type]] <- null_ss_df
}

## visualize
compare_res_ls <- list()
p_null <- list()
p_generank <- list()
colPtype <- c("HSC1-MPP1"="#0094D6", "HSC1-N2"="#F7AE24")
for (pair_type in c("HSC1-MPP1", "HSC1-N2")){
      real_curve <- ss_df_ls[[pair_type]]$ss_vector
      null_curves <- null_ss_df_ls[[pair_type]]
      compare_res <- compare_curve_no_align(real_curve = real_curve, null_curves = null_curves, 
                                          method = c("tss"), alternative = "more")
      compare_res_ls[[pair_type]] <- compare_res
      
      obs_score <- compare_res[["obs_score"]]
      null_scores <- compare_res[["null_scores"]]
      null_cutoff <- quantile(null_scores, probs = c(0.95))

      p_null[[pair_type]] <- ggplot(data.frame(x = null_scores), aes(x = x)) +
                                    geom_histogram(aes(y = ..density..), bins = 30, fill = "lightblue", color = "black", alpha = 0.6) +
                                    geom_density(color = "blue", size = 0.5) +
                                    geom_vline(xintercept = obs_score, color = "red", size = 0.5, linetype = 2) +
                                    geom_vline(xintercept = null_cutoff, color = "blue", size = 0.5, linetype = 2) +
                                    labs(title = "Null distribution", x = "The total sum of squares", y = "Density") +
                                    theme_classic()+
                                    theme(axis.text = element_text(size = 12, color = "black"), 
                                          axis.title = element_text(size = 12, color = "black"),
                                          plot.title = element_text(size = 12, color = "black"))+
                                    force_panelsizes(rows = unit(7, 'cm'), cols = unit(7, "cm"))

      ss_df <- ss_df_ls[[pair_type]]
      ss_df$Rank <- seq_along(ss_df$features)
      ss_df$Direction <- ifelse(ss_df$ss_vector > 0, "Up", "Dn")
      gene_label_df <- ss_df %>% dplyr::group_by(Direction) %>% dplyr::top_n(n = 10, wt = abs(ss_vector)) 
      p_generank[[pair_type]] <- ggplot(ss_df, aes(x = Rank, y = ss_vector)) +
                                    geom_point(size = 1, color = colPtype[pair_type])+
                                    geom_hline(yintercept = 0, linetype = 2, size = 0.5)+
                                    geom_text_repel(data = gene_label_df, aes(label=features), 
                                          color="black", size=3, 
                                          fontface="italic", arrow = arrow(ends="first", length = unit(0.01, "npc")), 
                                          box.padding = 0.1,  
                                          segment.color = 'black', segment.size = 0.2, force = 0.8, max.iter = 1e4)+
                                    labs(title = pair_type, x = "Gene rank", y = "Sum of pair log2FC") +
                                    theme_classic()+
                                    theme(axis.text = element_text(size = 12, color = "black"), 
                                          axis.title = element_text(size = 12, color = "black"),
                                          plot.title = element_text(size = 12, color = "black"), 
                                          legend.position = "none")+
                                    force_panelsizes(rows = unit(15, 'cm'), cols = unit(7, "cm"))+
                                    coord_cartesian(clip = "off")
}

compare_res_ls[["HSC1-MPP1"]][["p_value"]]
compare_res_ls[["HSC1-N2"]][["p_value"]]

pdf("plot/78.pdc.randomforest.iml-explained-125-hsc1mpp1n2-features.pair-log2fc-distribution.null-distribution.pdf", 
      width = 15, height = 18)
print(plot_grid(plotlist=c(p_null, p_generank), ncol = 2))
dev.off()

ss_df
paired_obj <- obj %>% subset(RF_predict_pairs == "HSC1-N2")
test <- my_compute_pairwise_log2FC_ss(seu = paired_obj, gene_list = c("Pf4", "Flt3"), 
                                                cell_pairs = unique(paired_obj$Pairs),   
                                                ident.1 = "HSC1", ident.2 = "N2")
test




### -----------------------------------------------------------------------------------------
### check FACS marker of predicted PDC-N2s
load("Rdata/02.pdc.allgrp.randomforest-predict.rdata")
pdc$RF_predict_labels <- gsub("MPP6", "N1", pdc$RF_predict_labels)
pdc$RF_predict_labels <- gsub("MPP5", "N2", pdc$RF_predict_labels)
pdc$RF_predict_labels <- gsub("MPP7", "LyI", pdc$RF_predict_labels)
table(pdc$Group,pdc$RF_predict_labels)
obj <- pdc %>% subset(RF_predict_labels %in% c("HSC1", "N2", "MPP1"))

# label_order <-  c("HSC1", "N1", "N2", "LyI", paste0("MPP", seq(4)))
label_order <-  c("HSC1", "N2", "MPP1")
obj$RF_predict_labels <- factor(obj$RF_predict_labels, levels = label_order)
table(obj$Group, obj$RF_predict_labels)

features <- c("Slamf1", "Cd48", "Cd34", "Flt3", "Procr")
DotPlot(obj, features = features, group.by = "RF_predict_labels", scale = FALSE) +
      scale_color_gradientn(colors = c("#222A8A", "#D3CAB1", "#B70335"))+
      scale_size(range = c(2, 6))+
      labs(x = "", y = "")+
      theme(axis.text.x = element_text(color = "black", face = "italic", angle = 45, hjust = 1),
            axis.text.y = element_text(color = "black"))+
      force_panelsizes(rows = unit(5, "cm"), cols = unit(5,"cm"))
# ggsave("plot/79.pdc.randomforest.FACS-markers.in-RF_predict_labels.pdf")
ggsave("plot/80.pdc.randomforest.FACS-markers.in-RF_predict_labels.pdf")

table(obj$RF_predict_labels, obj$Group)
DimPlot(obj, reduction = "umap.cca", group.by = "RF_predict_labels")


DotPlot(obj %>% subset(RF_predict_labels == "HSC1"), features = features, group.by = "Group") +
      scale_color_gradientn(colors = c("#222A8A", "#D3CAB1", "#B70335"))+
      scale_size(range = c(2, 6))+
      labs(x = "", y = "")+
      theme(axis.text.x = element_text(color = "black", face = "italic", angle = 45, hjust = 1),
            axis.text.y = element_text(color = "black"))+
      force_panelsizes(rows = unit(5, "cm"), cols = unit(5,"cm"))



load("/home/mayao/LABdata/MY/hsc1_hsc2/mpp/Rdata/RCombineHSC1.02.hsc-n1n2n3.combine-obj.rdata")
DotPlot(combine, features = features, group.by = "Group") +
      scale_color_gradientn(colors = c("#222A8A", "#D3CAB1", "#B70335"))+
      scale_size(range = c(2, 6))+
      labs(x = "", y = "")+
      theme(axis.text.x = element_text(color = "black", face = "italic", angle = 45, hjust = 1),
            axis.text.y = element_text(color = "black"))+
      force_panelsizes(rows = unit(5, "cm"), cols = unit(5,"cm"))


## module score
load("/home/mayao/genome_file/GeneSets/2017_NAR.2020_CSC.2021_Blood.gene_set_list.rdata")
obj <- my_AddModuleScore(obj, gene_set_list = gene_set_list, gsids = names(gene_set_list))
colnames(obj@meta.data)
boxdata <- FetchData(obj %>% subset(RF_predict_labels %in% c("HSC1", "MPP1", "N2")) , vars = c("Group_RF", "HSC", paste0("MPP", seq(5)) ))
plotdata <- tidyr::pivot_longer(boxdata, cols = -Group_RF, names_to = "Genesets", values_to = "Score")
ggplot(plotdata, aes(x = Group_RF, y = Score, fill = Group_RF))+
    geom_boxplot(color = "black", outlier.shape = NA, alpha = 1) + 
    facet_wrap( vars(Genesets), ncol = 3 )+
#     scale_fill_manual(values = colGroup)+
    theme_classic() +
    theme(legend.title = element_blank(), 
          axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("plot/80.pdc.randomforest.genesets.in-RF_predict_labels.boxplot.pdf")


### Check Flt3 and CD34
obj$Group_RF <- paste0(obj$Group, "-", obj$RF_predict_labels)

flt3_data <- FetchData(obj %>% subset(RF_predict_labels %in% c("HSC1", "MPP1", "N2")), vars = c("Flt3", "Group_RF"))
colnames(flt3_data)[2] <- "Group"
temp <- FetchData(combine %>% subset(Group  %in% c("HSC1", "N1", "N2", "Ly-I")), vars = c("Flt3", "Group"))
flt3_combine <- rbind(flt3_data,temp)
table(flt3_combine$Group)
flt3_combine$Group <- factor(flt3_combine$Group, levels = c("Fresh_HSC1-HSC1", "HSC1", "N1", "N2", "Ly-I", 
                              "ST-HSC1", "ST-MPP1", "ST-N2", 
                              "S12-HSC1", "S12-MPP1", "S12-N2"))

ggplot(flt3_combine, aes(Group, Flt3, color = Group))+
      geom_boxplot()+
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("plot/81.pdc.randomforest.Flt3.in-RF_predict_labels.boxplot.pdf")


cd34_data <- FetchData(obj %>% subset(RF_predict_labels %in% c("HSC1", "MPP1", "N2")), vars = c("Cd34", "Group_RF"))
colnames(cd34_data)[2] <- "Group"
temp <- FetchData(combine %>% subset(Group  %in% c("HSC1", "N1", "N2", "Ly-I")), vars = c("Cd34", "Group"))
cd34_combine <- rbind(cd34_data,temp)
cd34_combine$Group <- factor(cd34_combine$Group, levels = c("Fresh_HSC1-HSC1", "HSC1", "N1", "N2", "Ly-I", 
                              "ST-HSC1", "ST-MPP1", "ST-N2", 
                              "S12-HSC1", "S12-MPP1", "S12-N2"))

ggplot(cd34_combine, aes(Group, Cd34, color = Group))+
      geom_boxplot()+
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("plot/82.pdc.randomforest.Cd34.in-RF_predict_labels.boxplot.pdf")



### -----------------------------------------------------------------------------------------
### check specific markers between MPP1 and N2 pairs
load("Rdata/02.pdc.allgrp.randomforest-predict.rdata")
pdc$RF_predict_labels <- gsub("MPP6", "N1", pdc$RF_predict_labels)
pdc$RF_predict_labels <- gsub("MPP5", "N2", pdc$RF_predict_labels)
pdc$RF_predict_labels <- gsub("MPP7", "LyI", pdc$RF_predict_labels)
table(pdc$Group,pdc$RF_predict_labels)

label_order <-  c("HSC1", "N1", "N2", "LyI", paste0("MPP", seq(4)))
pdc$RF_predict_labels <- factor(pdc$RF_predict_labels, levels = label_order)
colRef <- setNames(colRef, nm = label_order)

## ------------------------------------------------------------------------
## Only pairs and HSC1, MPP1 and N2

obj <- pdc %>% subset(isPair=="Paired")
table(obj$Group, obj$RF_predict_labels)
Pair_type_ls <- c()
for (i in seq(1, nrow(obj@meta.data), 2)){
      Pair_type <- paste0(obj@meta.data[i, "RF_predict_labels"], "-", obj@meta.data[i+1, "RF_predict_labels"])
      Pair_type_ls <- c(Pair_type_ls, Pair_type)
}
obj$Pair_type <- rep(Pair_type_ls, each = 2)

obj$Pair_type <- gsub("HSC1-MPP1", "MPP1-HSC1", obj$Pair_type)
obj$Pair_type <- gsub("N2-MPP1", "MPP1-N2", obj$Pair_type)
obj$Pair_type <- gsub("HSC1-N2", "N2-HSC1", obj$Pair_type)
table(obj$Pair_type)

sub <- obj %>% subset(Pair_type %in% c("MPP1-HSC1", "MPP1-MPP1", "MPP1-N2", "N2-HSC1", "N2-N2"))
feature_data <- FetchData(sub, vars = c("Itga2b", "Itgb3", "Mpl", "RF_predict_labels", "Group", "Pairs", "Pair_type"))
head(feature_data)

ggplot(feature_data %>% subset(RF_predict_labels %in% c("MPP1")), aes(x = Pair_type, y = Itga2b, fill = RF_predict_labels))+
      geom_boxplot()+
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(feature_data %>% subset(RF_predict_labels %in% c("MPP1")), aes(x = Pair_type, y = Mpl, fill = RF_predict_labels))+
      geom_boxplot()+
      theme(axis.text.x = element_text(angle = 45, hjust = 1))


### Correlation
# cor_data <- FetchData(sub %>% subset(Pair_type == "MPP1-N2"), vars = c("RF_predict_labels", "Pairs", rownames(sub)))
# cor_data[1:5,1:5]
# cor_val_ls <- list()
# for (i in seq(3, ncol(cor_data)) ){
#       cor_1 <- cor_data[cor_data$RF_predict_labels == "MPP1", i]
#       cor_2 <- cor_data[cor_data$RF_predict_labels == "N2", i]
#       df <- data.frame(cor_1 = cor_1, cor_2 = cor_2)
#       df <- df[order(df$cor_1, decreasing = FALSE), ] 

#       gene_name <- colnames(cor_data)[i]
#       cor_val_ls[[gene_name]] <- cor(df$cor_1, df$cor_2, method = "spearman")
# }
# cor_val <- cor_val_ls %>% unlist()
# cor_val <- cor_val[!is.na(cor_val)]
# cor_val[cor_val > 0.7]
# cor_val["Cd53"]

# plot_data <- data.frame(MPP1 = cor_data[cor_data$RF_predict_labels == "MPP1", "Cd53"], 
#                         N2 = cor_data[cor_data$RF_predict_labels == "N2", "Cd53"], 
#                         Pairs = unique(cor_data$Pairs))
# plot_data <- plot_data[order(plot_data$MPP1, decreasing = FALSE), ] 
# plot_data_mut <- tidyr::pivot_longer(plot_data, cols = -Pairs, names_to = "Pair_ident", values_to = "Expression")
# plot_data_mut$Pairs <- factor(plot_data_mut$Pairs, levels = rev(unique(plot_data_mut$Pairs)))

# ggplot(plot_data_mut, aes(x = Pairs, y = Expression, color = Pair_ident))+
#       geom_point()+
#       theme(axis.text.x = element_text(angle = 45, hjust = 1))



