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


### Define color palettes
colGroup <- c("HSC1_Fresh" = "black",
              "HSC1_ST_48h" = "grey40", 
              "HSC1_S12_72h" = "grey80", 
              "EMPP6_Fresh" = "#FF80C0",
              "EMPP6_ST_48h" = "#C3AED2")
colCondition <- c("Fresh"="#918994", "Cultured"="#D6493E")
colCellpop <- c("HSC1"="black", "EMPP6"="#FF80C0")
colMedium <- c("Fresh"="#A1D99B","ST"="#2171B5", "S12"="#ED85B0")
colBatch <- c("Batch1"="#94BFE6" , "Batch2"="#EC979F" , "Batch3"="#61B768")

colRef <- c("HSC1"="black", 
            "EMPP6"="#FF80C0", "EMPP5"="#188386", "Ly-I"="#808080",
            "MPP1"="navajowhite3", "MPP2"="darkorange1","MPP3"="forestgreen","MPP4"="royalblue1"
              )

## set python options and load modules
python_path <- '/home/mayao/.conda/envs/seurat/bin/'
use_python(python_path)
py_run_string("import numpy as np")
py_run_string("import matplotlib.pyplot as pl")


### set path
setwd("/home/mayao/LABdata/MY/pdc_combine")

### set Seurat
options(Seurat.object.assay.version = "v5")


## *******************************************************************************************
## Bulid reference (Figure 7B)
## *******************************************************************************************
## Loading hsc data
load("/home/mayao/LABdata/MY/hsc1_hsc2/new/Rdata/02.hsc.seurat.clustered.markered.rdata")
count_hsc <- hsc@assays[["RNA"]]$counts
meta_hsc <- hsc@meta.data[,c("Group","Marker")] %>% as.data.frame()
meta_hsc$Name <- "HSC"
meta_hsc$CD135_median <- "Null"
meta_hsc$CD34_median <- "Null"
colnames(count_hsc) <- paste0("HSC.", colnames(count_hsc))
rownames(meta_hsc) <- paste0("HSC.", rownames(meta_hsc))
all(colnames(count_hsc) == rownames(meta_hsc))

## loading public data
load("/home/mayao/LABdata/DF/01.UMI.TPM.Meta.Homeostasis.Rdata") # count TPM meta.bg
dim(TPM) # 28579  1270
dim(count) # 28579  1270
dim(meta.bg) # 1270   11
table(meta.bg$Sorting.cell)
meta.bg$Marker <- meta.bg$Marker %>% as.vector()
meta.bg[meta.bg$Sorting.cell == "ST_HSC", ]$Marker %>% table()
# CD34-Flt3-Lin-Kit+Sca1+ CD34+Flt3+Lin-Kit+Sca1+ 
#                      22                      19 
meta.bg[meta.bg$Sorting.cell == "MPP1", ]$Marker %>% table()
# CD34-Flt3-Lin-Kit+Sca1+CD150+CD48- 
#                                72
meta.bg[meta.bg$Sorting.cell == "MPP2", ]$Marker %>% table()
# CD34-Flt3-Lin-Kit+Sca1+CD150+CD48+ 
#                                 72
meta.bg[meta.bg$Sorting.cell == "MPP3", ]$Marker %>% table()
# CD34-Flt3-Lin-Kit+Sca1+CD150-CD48+ 
#                                 76
meta.bg[meta.bg$Sorting.cell == "MPP4", ]$Marker %>% table()
# CD34+Flt3+Lin-Kit+Sca1+CD150-CD48+ 
#                                 78
meta.bg[meta.bg$Sorting.cell == "LT_HSC", ]$Marker %>% table()
# CD34-Flt3-Lin-Kit+Sca1+ 
#                      22

## Extract MPP1,2,3,4,ST_HSC,LT_HSC,LMPP
celltypes <- levels(meta.bg$Sorting.cell)
kepttypes <- celltypes[-which(celltypes %in% c("ST_HSC","LT_HSC","LMPP"))]
cells_df <- meta.bg %>% 
    subset(subset=Sorting.cell %in% kepttypes) %>% 
    rownames()
length(cells_df ) # 1176
count_df <- count[,cells_df ] %>% as.matrix()
dim(count_df) # 28579   1176
meta_df <- meta.bg[cells_df, ]
dim(meta_df) # 1176  11
meta_df <- meta_df[, c("Sorting.cell", "Marker")]
colnames(meta_df)[1] <- "Group"
meta_df$Name <- "DF"
meta_df$CD135_median <- "Null"
meta_df$CD34_median <- "Null"

## load MPPs
load("Rdata/02.mpp.seuratobj.rdata")
count_mpp <- mpp@assays[["RNA"]]$counts 
meta_mpp <- mpp@meta.data[,c("Group","Marker", "CD135_median", "CD34_median")] %>% as.data.frame()
meta_mpp$Name <- "MPP"
colnames(count_mpp) <- paste0("MPP.", colnames(count_mpp))
rownames(meta_mpp) <- paste0("MPP.", rownames(meta_mpp))
all(colnames(count_mpp) == rownames(meta_mpp))

## combine
count_ls <- c(list(count_df), list(count_hsc), list(count_mpp))
names(count_ls) <- c("DF", "HSC", "MPP")
meta <- rbind(meta_df, meta_hsc, meta_mpp)
meta$Group <- factor(meta$Group, levels=c("HSC1", "Fraction I", "Fraction II", "Fraction III",
                                          "ESLAM", "ESLAMSK",
                                          "HSC2","EMPP6", "EMPP5", "Ly-I",
                                          "MPP1","MPP2", "MPP3","MPP4",
                                          "HPC2", "HPC3",
                                          "CMP", "CLP","GMP", "MEP",
                                          "MK", "EryA", "EryB", 
                                          "B cell", "CD4T", "CD8T", "NK cell",
                                          "Granulocyte", "Monocyte", "Macrophage"))  

## Seurat obj
blood_new <- CreateSeuratObject(counts = count_ls, meta.data = meta, min.cells = 0, min.features = 0)
table(blood_new$Group)
save(blood_new, file="Rdata/temp13.blood_new.before_any_process.notfiltergenes.rdata")

## Process only HSC & MPP
load("Rdata/temp13.blood_new.before_any_process.notfiltergenes.rdata")
obj <- blood_new %>% subset(subset = Group %in% c("HSC1", "EMPP6", "EMPP5", "Ly-I", paste0("MPP", seq(4))))
obj$Group <- factor(obj$Group, levels = c("HSC1", "EMPP6", "EMPP5", "Ly-I", paste0("MPP", seq(4))))
table(obj$Group)

## pct_counts_Mito
obj <- PercentageFeatureSet(obj, 
			pattern = "^mt-", 
			col.name = "pct_counts_Mito")

## Decide parameters
assay <- "RNA"
normalization.method <- "LogNormalize"
nfeature <- 1500 
ndim <- 5 
integrate_method <- "jpca"
wt <- 25

## Normalize and CC
obj <- NormalizeData(obj, assay="RNA")
s.genes <- cc.genes$s.genes %>% str_to_title()
g2m.genes <- cc.genes$g2m.genes %>% str_to_title()
obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
obj <- CellCycleScoring(obj, g2m.features=g2m.genes, s.features=s.genes)
obj$cc_difference <- obj$S.Score - obj$G2M.Score
obj$Phase <- factor(obj$Phase, levels=c("G1","S","G2M"))
obj[[assay]] <- split(obj[[assay]], f=obj$Name)

## Integrate
obj <- FindVariableFeatures(obj, assay="RNA", nfeature = nfeature)
obj <- ScaleData(obj, assay="RNA", vars.to.regress=c('pct_counts_Mito',"cc_difference"))
obj <- RunPCA(obj, assay="RNA", npcs=30)
obj <- my_integrate(seuratobj = obj, assay="RNA", method=integrate_method, wt= wt, normalization.method = normalization.method, ndim=ndim)
obj[["RNA"]] <- JoinLayers(obj[["RNA"]])

## UMAP
reduction_name <- paste0("umap.", integrate_method)
obj <- RunUMAP(obj, 
            reduction = integrate_method, 
            reduction.name = reduction_name, 
            dims=1:ndim, 
            min.dist=0.5,
            n.neighbors=30,
            verbose=FALSE)

## Plot umap
colGroup <- c("HSC1"="black", "EMPP5"="#188386", "EMPP6"="#FF80C0", "Ly-I"="#808080",
              "MPP1"="#E7B471", "MPP2"="#674F9B", "MPP3"="#73955A", "MPP4"="#6885C1")
umap <- FetchData(obj, vars = c("umapjpca_1", "umapjpca_2", "seurat_clusters", "Group"))
ggplot(umap, aes(umapjpca_1, umapjpca_2, fill=Group))+
                    geom_point(size = 1.2, shape = 21, color = "white") +
                    scale_fill_manual(values = colGroup, name="") +
                    theme_few() +
                    labs(x="UMAP 1", y="UMAP 2") +
                    guides(color = guide_legend(override.aes = list(size=3))) +
                    force_panelsizes(rows=unit(6,"cm"), cols=unit(6,"cm"))
ggsave(paste0("plot/17-3.hspc-combine.260206_new.nfeature-",nfeature,".",integrate_method,"-",ndim,"-",wt,".umap.pdf"))

## save rdata
mu_hspc <- obj
saveRDS(mu_hspc, file = "Rdata/03.combine.hspc.seuratobj.rdata")



## *******************************************************************************************
## Train model (Figure 7C)
## *******************************************************************************************
load("Rdata/03.combine.hspc.seuratobj.rdata")
obj <- hspc

## random forest
library(caret)
library(randomForest)
rf_df <- FetchData(obj, vars = c(VariableFeatures(obj), "Group")) # 1500
dim(rf_df)
colnames(rf_df) <- gsub("-", "_", colnames(rf_df))
colnames(rf_df) <- paste0("G_", colnames(rf_df) ) # so that no colnames are begining with number
colnames(rf_df)[ncol(rf_df)] <- "Group"

train_index <- createDataPartition(rf_df$Group, p = 0.7, list = FALSE)
train_data <- rf_df[train_index, ]
test_data <- rf_df[-train_index, ]
rf_model_A <- randomForest(Group ~ ., data = train_data, importance = TRUE)
predictions <- predict(rf_model_A, test_data)
result_df <- data.frame(row.names = rownames(test_data), 
                  True_labels = test_data$Group, 
                  Predict_labels = predictions)

## importance of featuress
score_df <- importance(rf_model_A) %>% as.data.frame()
score_df$Score <- score_df$MeanDecreaseAccuracy + score_df$MeanDecreaseGini
score_df <- score_df[order(score_df$Score, decreasing = TRUE), ]
score_df$E5_E6 <- score_df$EMPP5 * score_df$EMPP6
score_df$E5_E6_LyI <- score_df$EMPP5 * score_df$EMPP6 * score_df$Ly-I
quantile(score_df$E5_E6_LyI, seq(0, 1, 0.1))
selected <- rownames(score_df)[score_df$E5_E6 > 0.15155351] 
length(selected) # 300


## re-model
train_data_B <- subset(train_data, Group %in% c("EMPP5", "EMPP6", "Ly-I"))
train_data_B$Group <- as.character(train_data_B$Group) %>% factor()
table(train_data_B$Group)
rf_model_B <- randomForest(Group ~ ., data = train_data_B[, c(selected, "Group")], importance = TRUE)

table(result_df$Predict_labels)
test_data_cells_B <- result_df %>% subset(Predict_labels %in% c("EMPP5", "EMPP6", "Ly-I")) %>% rownames()
test_data_B<- test_data[test_data_cells_B, c(selected, "Group")]
test_data_B$Group <- as.character(test_data_B$Group) %>% factor()
table(test_data_B$Group)
predictions_B <- predict(rf_model_B , test_data_B)
result_df_B <- data.frame(row.names = rownames(test_data_B), 
                  True_labels = test_data_B$Group, 
                  Predict_labels = predictions_B)
table(result_df_B$True_labels, result_df_B$Predict_labels)

result_df$Predict_labels_B <- sapply(rownames(result_df), function(x){
    if (result_df[x, "True_labels"] %in% c("EMPP5", "EMPP6", "Ly-I")){
        Predict_labels_B <- result_df_B[x, "Predict_labels"]
    }else {
        Predict_labels_B <- result_df[x, "Predict_labels"]
    }
    return(Predict_labels_B)
})

table(result_df$True_labels, result_df$Predict_labels_B)
backup <- result_df
result_df$Predict_labels <- result_df$Predict_labels_B
result_df$Predict_labels_B <- NULL
table(result_df$True_labels, result_df$Predict_labels)


accuracy <- as.vector(table(result_df$True_labels == result_df$Predict_labels)/nrow(result_df))[2] %>% 
                    signif(3)

test_df_c <- ddply(result_df, .(True_labels, Predict_labels), summarise, Count = length(Predict_labels))
test_df_p <- ddply(test_df_c, .(True_labels ), summarise, Predict_labels=Predict_labels, Percent = Count/sum(Count)*100)
ggplot(test_df_p, aes(x = factor(Predict_labels), y = factor(True_labels), fill = Percent)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "#E41A1C") +
    labs(title = paste0(""), x = "Predicted Label", y = "Original Label") +
    theme_minimal()+
    theme(axis.text.x = element_text(angle=45,hjust=1))+
    ggtitle("")+
    force_panelsizes(rows=unit(7,"cm"), cols=unit(7,"cm"))
ggsave("plot/23.hspc-combine.randomforest.test-tile-plot.pdf")


### Save model
features <- VariableFeatures(obj)
selected_features_for_B <- selected
save(rf_model_A, rf_model_B, features, selected_features_for_B, file = "Rdata/04.hspc-combine.rf_model.randomforest.rdata")
saveRDS(rf_model_A, "Rdata/04.hspc-combine.rf_model_A.randomforest.rds")
saveRDS(rf_model_B, "Rdata/04.hspc-combine.rf_model_B.randomforest.rds")
saveRDS(features, "Rdata/04.hspc-combine.features.randomforest.rds")
saveRDS(selected_features_for_B, "Rdata/04.hspc-combine.selected_features_for_B.randomforest.rds")



## *******************************************************************************************
## Predict HSC1-PDC (Figure 7E, S7C, S7E)
## *******************************************************************************************
load("/home/mayao/LABdata/MY/pdc/Rdata/temp03.pdc.seurat.dim_reduc.clustered.rdata")
load("Rdata/04.hspc-combine.rf_model.randomforest.rdata")
query <- pdc %>% subset(subset = Name != "Exp5")
dim(query)
table(query$Group)

temp <- matrix(0, nrow = length(features[!features %in% rownames(query)]), ncol = ncol(query))
rownames(temp) <- features[!features %in% rownames(query)]
colnames(temp) <- colnames(query)

query_data <- rbind(query[["RNA"]]$data[features[features %in% rownames(query)],], temp) %>% t() %>% as.data.frame()
dim(query_data) # 640 1500

query_data$Group <- query@meta.data[match(rownames(query_data), rownames(query@meta.data)), "Group"]
dim(query_data) # 640 1501

colnames(query_data) <- gsub("-", "_", colnames(query_data))
colnames(query_data) <- paste0("G_", colnames(query_data) ) # so that no colnames are begining with number
colnames(query_data)[ncol(query_data)] <- "Group"
table(query_data$Group)

predictions_A <- predict(rf_model_A, query_data)
predictions_A_score <- predict(rf_model_A, query_data, type = "prob")
res_A <- data.frame(row.names = rownames(query_data), 
                  Group = query_data$Group, 
                  Predict_labels = predictions_A)
pheatmap::pheatmap(predictions_A_score)
predictions_A_score[,"MPP1"] %>% max()

query_data_cells_B <- res_A %>% subset(Predict_labels %in% c("EMPP5", "EMPP6", "Ly-I")) %>% rownames()
query_data_B<- query_data[query_data_cells_B, c(selected_features_for_B, "Group")]
table(query_data_B$Group)
predictions_B <- predict(rf_model_B , query_data_B)
predictions_B_score <- predict(rf_model_B, query_data_B, type = "prob")
res_B <- data.frame(row.names = rownames(query_data_B), 
                  True_labels = query_data_B$Group, 
                  Predict_labels = predictions_B)
table(res_B$True_labels, res_B$Predict_labels)
pheatmap::pheatmap(predictions_B_score)


pdc_res <- data.frame(row.names = rownames(query_data), 
                  Group = query_data$Group)

pdc_res$Predict_labels <- sapply(rownames(res_A), function(x){
    if (res_A[x, "Predict_labels"] %in% c("EMPP5", "EMPP6", "Ly-I")){
        Predict_labels <- res_B[x, "Predict_labels"]
    }else {
        Predict_labels <- res_A[x, "Predict_labels"]
    }
    return(Predict_labels)
})
predict_label_order <- c("HSC1", "EMPP6", "EMPP5", "Ly-I", paste0("MPP", seq(4)))
pdc_res$Predict_labels <- factor(pdc_res$Predict_labels, levels = predict_label_order)

save(pdc, file = "/home/mayao/LABdata/MY/pdc/Rdata/02.pdc.allgrp.randomforest-predict.rdata")



## *******************************************************************************************
## Predict EMPP6-PDC (Figure 7E)
## *******************************************************************************************
load("Rdata/04.pdc_E6.seurat.dim_reduc.clustered.rdata")
load("/home/mayao/LABdata/MY/hsc1_hsc2/mpp/Rdata/04.hspc-combine.rf_model.randomforest.rdata")
query <- pdc_E6
dim(query)
table(query$Condition)

temp <- matrix(0, nrow = length(features[!features %in% rownames(query)]), ncol = ncol(query))
rownames(temp) <- features[!features %in% rownames(query)]
colnames(temp) <- colnames(query)

query_df <- rbind(query[["RNA"]]$data[features[features %in% rownames(query)],], temp) %>% t() %>% as.data.frame()
dim(query_df) # 359 1500

query_df$Group <- query@meta.data[match(rownames(query_df), rownames(query@meta.data)), "Condition"]
dim(query_df) # 359 1501

colnames(query_df) <- gsub("-", "_", colnames(query_df))
colnames(query_df) <- paste0("G_", colnames(query_df) ) # so that no colnames are begining with number
colnames(query_df)[ncol(query_df)] <- "Group"
table(query_df$Group)

query_data <- subset(query_df, select = -Group)
predictions_A <- predict(rf_model_A, query_data)
res_A <- data.frame(row.names = rownames(query_df), 
                  Group = query_df$Group, 
                  Predict_labels = predictions_A)

query_data_cells_B <- res_A %>% subset(Predict_labels %in% c("EMPP5", "EMPP6", "Ly-I")) %>% rownames()
query_data_B<- query_data[query_data_cells_B, c(selected_features_for_B)]
predictions_B <- predict(rf_model_B , query_data_B)
res_B <- data.frame(row.names = rownames(query_data_B), 
                  True_labels = query_df[query_data_cells_B,]$Group, 
                  Predict_labels = predictions_B)


pdc_res <- data.frame(row.names = rownames(query_df), 
                  Group = query_df$Group)

pdc_res$Predict_labels <- sapply(rownames(res_A), function(x){
    if (res_A[x, "Predict_labels"] %in% c("EMPP5", "EMPP6", "Ly-I")){
        Predict_labels <- res_B[x, "Predict_labels"]
    }else {
        Predict_labels <- res_A[x, "Predict_labels"]
    }
    return(Predict_labels)
})
pdc_res$Predict_labels <- factor(pdc_res$Predict_labels, levels = predict_label_order)
pdc_E6$RF_predict_labels <-  pdc_res$Predict_labels[match(rownames(pdc_res), rownames(pdc_E6@meta.data))]
save(pdc_E6, file = "/home/mayao/LABdata/MY/pdc_E6/Rdata/05.pdc_E6.randomforest-predict.rdata")





## *******************************************************************************************
## Combine HSC1-PDC and EMPP6-PDC
## Note: Whether HSC-PDC and E6-PDC were predicted separately or together,
##       RF predict results remained the same.
## *******************************************************************************************
load("/home/mayao/LABdata/MY/pdc/Rdata/02.pdc.allgrp.randomforest-predict.rdata")
load("/home/mayao/LABdata/MY/pdc_E6/Rdata/05.pdc_E6.randomforest-predict.rdata")

## Add bacth info
pdc$Batch <- "Batch1"
pdc_E6$Batch <- "Batch2"
colnames(pdc_E6@meta.data)

## HSC1 pdc meta
pdc_meta <- data.frame(row.names = rownames(pdc@meta.data), Group = pdc$Group, Condition = pdc$Condition, Cell_pop = "HSC1", Pairs = pdc$Pairs, Name = pdc$Batch)
pdc_meta$Group <- gsub("Fresh_HSC1", 'HSC1_Fresh', pdc_meta$Group, fixed = TRUE)
pdc_meta$Group <- gsub("ST", 'HSC1_ST_48h', pdc_meta$Group, fixed = TRUE)
pdc_meta$Group <- gsub("S12", 'HSC1_S12_72h', pdc_meta$Group, fixed = TRUE)
table(pdc_meta$Group)

## N1 pdc meta
pdc_E6_meta <- data.frame(row.names = rownames(pdc_E6@meta.data), Group = pdc_E6$Group, Condition = pdc_E6$Condition, Cell_pop = "N1", Pairs = pdc_E6$Pairs, Name = pdc_E6$Batch)
pdc_E6_meta$Group <- gsub("Fresh", 'N1_Fresh', pdc_E6$Condition)
pdc_E6_meta$Group <- gsub("ST", 'N1_ST_48h', pdc_E6_meta$Group)
table(pdc_E6_meta$Group)

## get integrated meta
meta_ls <- c(list(pdc_meta), list(pdc_E6_meta))
meta <- Reduce(rbind, meta_ls)
meta$Pairs <- gsub("Fresh", "Not_pairs", meta$Pairs)
meta$Condition <- gsub("ST", "Cultured", meta$Condition)
meta$Condition <- gsub("S12", "Cultured", meta$Condition)
meta$Culture_time <- str_split(meta$Group, "[_]", simplify = TRUE)[,3]
meta$Culture_time[meta$Culture_time == ""] <- "0h"
meta$Medium <- str_split(meta$Group, "[_]", simplify = TRUE)[,2]
table(meta$Condition)
table(meta$Group)
table(meta$Pairs)
table(meta$Cell_pop)
table(meta$Name)
table(meta$Culture_time)

## get count_ls
count_ls <- c(list(pdc[["RNA"]]$counts), list(pdc_E6[["RNA"]]$counts))

## Seurat
obj <- CreateSeuratObject(counts = count_ls, meta.data = meta, min.cells = 0, min.features = 0)

## Decide parameters
assay <- "RNA"
normalization.method <- "LogNormalize"
nfeature <- 2500
ndim <- 6
integrate_method <- "cca"
wt <- 10

## log norm
obj <- NormalizeData(obj, assay="RNA")

## cc genes
s.genes <- cc.genes$s.genes %>% str_to_title()
g2m.genes <- cc.genes$g2m.genes %>% str_to_title()

## join layers
obj[["RNA"]] <- JoinLayers(obj[["RNA"]])

## Mito and CC
obj <- PercentageFeatureSet(obj, 
			pattern = "^mt-", 
			col.name = "pct_counts_Mito")
obj <- CellCycleScoring(obj, g2m.features=g2m.genes, s.features=s.genes)
obj$cc_difference <- obj$S.Score - obj$G2M.Score
obj$Phase <- factor(obj$Phase, levels=c("G1","S","G2M"))
table(obj$Phase, obj$Group)

## Split
obj[[assay]] <- split(obj[[assay]], f=obj$Name)

## Variable features
obj <- FindVariableFeatures(obj, assay="RNA", nfeature = nfeature)

## scale
obj <- ScaleData(obj, assay="RNA", vars.to.regress=c('pct_counts_Mito'))

## PCA
obj <- RunPCA(obj, assay="RNA", npcs=30)

# Integrate
obj <- my_integrate(seuratobj = obj, assay="RNA", method=integrate_method, wt= wt, normalization.method = normalization.method, ndim=ndim)
obj[["RNA"]] <- JoinLayers(obj[["RNA"]])

## UMAP
reduction_name <- paste0("umap")
obj <- RunUMAP(obj, 
            reduction = integrate_method, 
            reduction.name = reduction_name, 
            dims=1:ndim, 
            verbose=FALSE)
pdc_combine <- obj 
save(pdc_combine, file = "Rdata/03.pdc_nocontrol.seuratobj.rdata")

data_df <- FetchData(obj, vars = c("umaprpca_1", "umaprpca_2", "Cell_pop", "Culture_time", "Medium", "seurat_clusters", "Group"))
data_df$Cell <- rownames(data_df)
head(data_df)


## ********** Figure 7D **********
colGroup <- c("HSC1_Fresh" = "black",
              "HSC1_ST_48h" = "grey40", 
              "HSC1_S12_72h" = "grey80", 
              "N1_Fresh" = "#FF80C0",
              "N1_ST_48h" = "#C3AED2")
ggplot(data_df, aes(umaprpca_1, umaprpca_2, color = Group))+
  geom_point(size = 0.8)+
  scale_color_manual(values = colGroup)+
  theme_classic()+
  theme(axis.text = element_text(color = "black", size = 10))+
  guides(colour = guide_legend(override.aes = list(size=3), title=NULL) ) +
  labs(title = "", x = "UMAP 1", y = "UMAP 2")+
  force_panelsizes(rows = unit(7,"cm"), cols = unit(7,"cm"))
ggsave("plot/33.pdc_nocontrol.umap.pdf", width = 10, height = 14)
