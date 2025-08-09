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

## predict label order
predict_label_order <- c("HSC1", paste0("MPP", seq(7)))


## --------------------------------------------------------------------------------------------
load("Rdata/03.combine.hspc.seuratobj.rdata")
obj <- hspc

## ---------------------------------------------------------------------------------------------
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
table(result_df$True_labels, result_df$Predict_labels)


## importance of featuress
score_df <- importance(rf_model_A) %>% as.data.frame()
score_df$Score <- score_df$MeanDecreaseAccuracy + score_df$MeanDecreaseGini
score_df <- score_df[order(score_df$Score, decreasing = TRUE), ]
score_df$MPP5_MPP6 <- score_df$MPP5 * score_df$MPP6
score_df$MPP5_MPP6_MPP7 <- score_df$MPP5 * score_df$MPP6 * score_df$MPP7
quantile(score_df$MPP5_MPP6_MPP7, seq(0, 1, 0.1))
selected <- rownames(score_df)[score_df$MPP5_MPP6 > 0.15155351] 
length(selected) # 300


## re-model
train_data_B <- subset(train_data, Group %in% c("MPP5", "MPP6", "MPP7"))
train_data_B$Group <- as.character(train_data_B$Group) %>% factor()
table(train_data_B$Group)
rf_model_B <- randomForest(Group ~ ., data = train_data_B[, c(selected, "Group")], importance = TRUE)

table(result_df$Predict_labels)
test_data_cells_B <- result_df %>% subset(Predict_labels %in% c("MPP5", "MPP6", "MPP7")) %>% rownames()
test_data_B<- test_data[test_data_cells_B, c(selected, "Group")]
test_data_B$Group <- as.character(test_data_B$Group) %>% factor()
table(test_data_B$Group)
predictions_B <- predict(rf_model_B , test_data_B)
result_df_B <- data.frame(row.names = rownames(test_data_B), 
                  True_labels = test_data_B$Group, 
                  Predict_labels = predictions_B)
table(result_df_B$True_labels, result_df_B$Predict_labels)

result_df$Predict_labels_B <- sapply(rownames(result_df), function(x){
    if (result_df[x, "True_labels"] %in% c("MPP5", "MPP6", "MPP7")){
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

## calculate TP, FP, TN, FN
temp_c <- ddply(result_df, .(True_labels, Predict_labels),summarize, Count=length(Predict_labels))
temp_c %>% head(n = 10)
temp_p <- ddply(temp_c, .(True_labels),summarize, Predict_labels=Predict_labels, Percent=Count*100/sum(Count))
temp_p %>% head(n = 10)
temp_p$Percent <- signif(temp_p$Percent, 3)

true_pos <- temp_p[temp_p$True_labels == temp_p$Predict_labels,]
false_neg <- temp_p[temp_p$True_labels != temp_p$Predict_labels,] %>% 
                plotly::group_by(True_labels) %>% 
                plotly::summarise(True_labels=True_labels, FN = sum(Percent)) %>% 
                distinct()
p_fn <- ggplot(temp_p, aes(x = factor(Predict_labels), y = factor(True_labels), fill = Percent)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "blue") +
    labs(title = paste0(""), x = "Predicted Label", y = "Original Label") +
    theme_minimal()+
    theme(axis.text.x = element_text(angle=45,hjust=1))+
    ggtitle("False negative")+
    geom_text(aes(label=Percent))+
    force_panelsizes(rows=unit(7,"cm"), cols=unit(7,"cm"))
p_fn

temp_c <- ddply(result_df, .(Predict_labels, True_labels),summarize, Count=length(True_labels))
temp_p <- ddply(temp_c, .(Predict_labels),summarize, True_labels=True_labels, Percent=Count*100/sum(Count))
temp_p$Percent <- signif(temp_p$Percent, 1)
false_pos <- temp_p[temp_p$Predict_labels != temp_p$True_labels,] %>% 
                plotly::group_by(Predict_labels) %>% 
                plotly::summarise(Predict_labels=Predict_labels, FP = sum(Percent)) %>% 
                distinct()
p_fp <- ggplot(temp_p, aes(x = factor(Predict_labels), y = factor(True_labels), fill = Percent)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "red") +
    labs(title = paste0(""), x = "Predicted Label", y = "Original Label") +
    theme_minimal()+
    theme(axis.text.x = element_text(angle=45,hjust=1))+
    ggtitle("False positive")+
    geom_text(aes(label=Percent))+
    force_panelsizes(rows=unit(7,"cm"), cols=unit(7,"cm"))
p_fp

accuracy <- as.vector(table(result_df$True_labels == result_df$Predict_labels)/nrow(result_df))[2] %>% 
                    signif(3)
accuracy 

test_df_c <- ddply(result_df, .(True_labels, Predict_labels), summarise, Count = length(Predict_labels))
test_df_p <- ddply(test_df_c, .(True_labels ), summarise, Predict_labels=Predict_labels, Percent = Count/sum(Count)*100)
ggplot(test_df_p, aes(x = factor(Predict_labels), y = factor(True_labels), fill = Percent)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "#E41A1C") +
    labs(title = paste0(""), x = "Predicted Label", y = "Original Label") +
    theme_minimal()+
    theme(axis.text.x = element_text(angle=45,hjust=1))+
    ggtitle("")+
    # geom_text(aes(label=Percent))+
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



###### ----------------------------------------------------------------------------------------
###### Predict all HSPC data
library(randomForest)
packageVersion("randomForest")
load("Rdata/04.hspc-combine.rf_model.randomforest.rdata")
## --------------------------------------------------------------------------------------------
load("Rdata/03.combine.hspc.seuratobj.rdata")
query <- hspc
dim(query)
table(query$Group)

# query_data <- rbind(query[["RNA"]]$data[features,], temp) %>% t() %>% as.data.frame()
query_data <- query[["RNA"]]$data[features,] %>% t() %>% as.data.frame()
dim(query_data) # 696 1500

query_data$Group <- query@meta.data[match(rownames(query_data), rownames(query@meta.data)), "Group"]
dim(query_data) # 696 1501

colnames(query_data) <- gsub("-", "_", colnames(query_data))
colnames(query_data) <- paste0("G_", colnames(query_data) ) # so that no colnames are begining with number
colnames(query_data)[ncol(query_data)] <- "Group"
table(query_data$Group)

predictions_A <- predict(rf_model_A, query_data)
res_A <- data.frame(row.names = rownames(query_data), 
                  Group = query_data$Group, 
                  Predict_labels = predictions_A)

query_data_cells_B <- res_A %>% subset(Predict_labels %in% c("MPP5", "MPP6", "MPP7")) %>% rownames()
query_data_B<- query_data[query_data_cells_B, c(selected_features_for_B, "Group")]
table(query_data_B$Group)
predictions_B <- predict(rf_model_B , query_data_B)
res_B <- data.frame(row.names = rownames(query_data_B), 
                  True_labels = query_data_B$Group, 
                  Predict_labels = predictions_B)
table(res_B$True_labels, res_B$Predict_labels)

predict_res <- data.frame(row.names = rownames(query_data), 
                  Group = query_data$Group)

predict_res$Predict_labels <- sapply(rownames(res_A), function(x){
    if (res_A[x, "Predict_labels"] %in% c("MPP5", "MPP6", "MPP7")){
        Predict_labels <- res_B[x, "Predict_labels"]
    }else {
        Predict_labels <- res_A[x, "Predict_labels"]
    }
    return(Predict_labels)
})
table(predict_res$Group, predict_res$Predict_labels)
predict_res$Predict_labels <- factor(predict_res$Predict_labels, levels = predict_label_order)
table(predict_res$Group, predict_res$Predict_labels)
nrow(predict_res)

## change data frame
result_df <- predict_res
colnames(result_df) <- c("True_labels", "Predict_labels")

## confusion matrix
library(ggplot2)
library(reshape2)
library(scales)
result_df$True_labels <- gsub("MPP5", "N2", result_df$True_labels)
result_df$True_labels <- gsub("MPP6", "N1", result_df$True_labels)
result_df$True_labels <- gsub("MPP7", "N3", result_df$True_labels)
result_df$True_labels <- factor(result_df$True_labels, levels = c("HSC1", "N1", "N2", "N3", paste0("MPP", seq(4))))

result_df$Predict_labels <- gsub("MPP5", "N2", result_df$Predict_labels)
result_df$Predict_labels <- gsub("MPP6", "N1", result_df$Predict_labels)
result_df$Predict_labels <- gsub("MPP7", "N3", result_df$Predict_labels)
result_df$Predict_labels <- factor(result_df$Predict_labels, levels = c("HSC1", "N1", "N2", "N3", paste0("MPP", seq(4))))

conf_matrix <- table(result_df$True_labels, result_df$Predict_labels)
conf_matrix_percent <- prop.table(conf_matrix, 1) * 100
TPR <- diag(conf_matrix_percent)
FNR <- 100 - TPR
df_melt <- melt(conf_matrix_percent)
tpr_fnr_df <- data.frame(
  True_labels = rownames(conf_matrix_percent),
  TPR = TPR,
  FNR = FNR
)
ggplot(df_melt, aes(x = Var2, y = Var1)) +
    geom_tile(aes(fill = value), color = "white") +
    scale_fill_gradientn(colors = colorRampPalette(c("white","#82CCD5", "#EF8FA2"))(100), limits = c(0, 100))+
    #   scale_fill_gradient(low = "white", high = "green", limits = c(0, 100)) +
    geom_text(aes(label = paste0(round(value, 1), "%")), color = "black", size = 3) +
    theme_bw() +
    labs(x = "Predicted class", y = "True class", fill = "Percent") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), 
          panel.grid.major = element_blank()) +
    ggtitle("Accuracy = 96.1%")+
    force_panelsizes(rows=unit(6.5,"cm"), cols=unit(6.5,"cm"))
ggsave("plot/23.hspc-combine.randomforest.test-tile-plot.pdf")




###### ----------------------------------------------------------------------------------------
###### PDC Predictions
library(randomForest)
packageVersion("randomForest")
load("/home/mayao/LABdata/MY/pdc/Rdata/temp03.pdc.seurat.dim_reduc.clustered.rdata")
load("Rdata/04.hspc-combine.rf_model.randomforest.rdata")
query <- pdc %>% subset(subset = Name != "Exp5")
dim(query)
table(query$Group)

# query[["RNA"]]$scale.data <- NULL
# query <- ScaleData(query, features = features, assay = "RNA")

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

query_data_cells_B <- res_A %>% subset(Predict_labels %in% c("MPP5", "MPP6", "MPP7")) %>% rownames()
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
    if (res_A[x, "Predict_labels"] %in% c("MPP5", "MPP6", "MPP7")){
        Predict_labels <- res_B[x, "Predict_labels"]
    }else {
        Predict_labels <- res_A[x, "Predict_labels"]
    }
    return(Predict_labels)
})
table(pdc_res$Group, pdc_res$Predict_labels)
pdc_res$Predict_labels <- factor(pdc_res$Predict_labels, levels = predict_label_order)
table(pdc_res$Group, pdc_res$Predict_labels)
nrow(pdc_res)


df_c <- ddply(pdc_res, .(Group, Predict_labels), summarise, Count = length(Predict_labels))
df_p <- ddply(df_c, .(Group), summarise, Predict_labels=Predict_labels, Percent = Count/sum(Count)*100)

ggplot(df_p, aes(x = Group, y = Percent, fill = Predict_labels, stratum = Predict_labels, alluvium =  Predict_labels)) + 
      geom_flow(width = 0.5, alpha = 0.3, knot.pos=0, color = 'white') + 
      geom_col(width = 0.5, color = 'white') + 
      scale_y_continuous(expand = c(0, 0)) +  
      scale_fill_manual(values = colRef, name = "") +
      xlab("") + 
      ylab("% Pairs") + 
      theme_classic() +
      theme(legend.title = element_blank())+
      force_panelsizes(rows=unit(6,"cm"), cols=unit(5,"cm"))
pdc <- query
dim(pdc)
pdc$RF_predict_labels <-  pdc_res$Predict_labels[match(rownames(pdc_res), rownames(pdc@meta.data))]
# save(pdc, file = "/home/mayao/LABdata/MY/pdc/Rdata/02.pdc.allgrp.randomforest-predict.rdata")



###### ----------------------------------------------------------------------------------------
library(caret)
library(randomForest)

## load data
load("Rdata/04.hspc-combine.rf_model.randomforest.rdata")
load("Rdata/03.combine.hspc.seuratobj.rdata")
obj <- hspc
rf_df <- FetchData(obj, vars = c(features, "Group")) # 1500
colnames(rf_df) <- gsub("-", "_", colnames(rf_df))
colnames(rf_df) <- paste0("G_", colnames(rf_df) ) # so that no colnames are begining with number
colnames(rf_df)[ncol(rf_df)] <- "Group"

### ----------------------------------------------------------------------------------------------
### iml explain
library(iml)
library(mlr)


## 1. Get top genes with highest influence (but without direction)
get_top_genes_for_all_classes <- function(rf_model, X, y, top_n = 10) {
  if (!requireNamespace("iml", quietly = TRUE)) install.packages("iml")
  library(iml)

  X <- as.data.frame(X)
  y <- as.factor(y)

  # get prob matrix of all classes
  prob_matrix <- predict(rf_model, X, type = "prob")
  all_classes <- colnames(prob_matrix)

  imp_list <- list()
  gene_list <- list()

  for (class_name in all_classes) {
    message(paste0("Calculating feature importance for [", class_name, "]..."))

    # specific class
    class_prob <- prob_matrix[, class_name]

    # predict function for iml
    predict_class_prob <- function(model, newdata) {
      newdata <- as.data.frame(newdata)
      predict(model, newdata = newdata, type = "prob")[, class_name]
    }

    # predictor
    predictor <- Predictor$new(
      model = rf_model,
      data = X,
      y = class_prob,
      predict.function = predict_class_prob,
      type = "response"
    )

    # importance
    imp <- FeatureImp$new(predictor, loss = "mse")
    imp_list[[class_name]] <- imp

    # top genes
    top_genes <- head(imp$results[order(imp$results$importance, decreasing = TRUE), ], top_n)
    gene_list[[class_name]] <- top_genes
  }

  result_list <- c(list(imp_list), list(gene_list))
  names(result_list) <- c("imp_list", "gene_list")

  return(result_list)
}


res_ls <- get_top_genes_for_all_classes(
  rf_model = rf_model_A, 
  X = rf_df[, -ncol(rf_df)],
  y = rf_df$Group,
  top_n = 50
)

save(res_ls, file = "Rdata/09.1.hspc-combine.rf_model_A.iml-explanation.res_ls.rdata")



load("Rdata/09.1.hspc-combine.rf_model_A.iml-explanation.res_ls.rdata")
temp <- c(list(res_ls[1:8]), list(res_ls[9:16]))
names(temp) <- c("imp_list", "gene_list")

temp[["gene_list"]][["MPP6"]]

xlsx_file <- paste0("22.hspc.random-forest.iml.important-features.xlsx")
if ( file.exists(paste0("table/",xlsx_file)) ) {
  file.remove( paste0("table/",xlsx_file) )
}
wb <- my_xlsx_create(xlsx_file)

for (top_n in c(10, 30, 50, 100)){

  top_genes_ls <- list()
  for (cell_type in names(temp[["imp_list"]])){
    res_df <- temp[["imp_list"]][[cell_type]]$results
    res_df <- res_df[order(res_df$importance, decreasing = TRUE),]
    top_genes <- gsub("G_", "", res_df$feature)[1:top_n]
    top_genes <- gsub("_", "-", top_genes)
    top_genes_ls[[cell_type]] <- c(cell_type, top_genes)
  }
  top_genes_df <- Reduce(cbind, top_genes_ls)
  wb <- my_xlsx_write(wb, data = top_genes_df, sheet = paste0("top_",top_n), rowNames = TRUE)

}

my_xlsx_save(wb, xlsx_file)



### 2. Get top genes with highest influence with direction
get_positive_contributors_for_class <- function(rf_model, class_name, data, feature_to_fit) {
  library(iml)
  y_class <- predict(rf_model, data, type = "prob")[, class_name]
  predict_class_prob <- function(model, newdata) {
    predict(model, newdata = newdata, type = "prob")[, class_name]
  }
  predictor <- Predictor$new(
    model = rf_model,
    data = data,
    y = y_class,
    predict.function = predict_class_prob,
    type = "response"
  )

  # get slopes
  slopes <- numeric(length(feature_to_fit))

  message("Calculating slopes for [", cell_type, "]...")
  for (i in seq_along(feature_to_fit)) {
    message(" # Gene No.", i)
    fe <- FeatureEffect$new(predictor, feature = feature_to_fit[i], method = "pdp")

    eff_data <- fe$results
    colnames(eff_data)[1] <- "feature"
    fit <- lm(.value ~ feature, data = eff_data)
    slopes[i] <- coef(fit)[2]  
  }

  result <- data.frame(feature = feature_to_fit, slope = slopes)
  result <- result[order(result$slope, decreasing = TRUE), ]
  return(result)
}


load("Rdata/09.hspc-combine.rf_model_A.iml-explanation.res_ls.rdata")
temp <- c(list(res_ls[1:8]), list(res_ls[9:16]))
names(temp) <- c("imp_list", "gene_list")
plot( temp[["imp_list"]][["HSC1"]])


top_n <- 100
rf_data <- rf_df[, -ncol(rf_df)] %>% as.data.frame()
genes_res_ls <- list()
for (cell_type in names(temp[["imp_list"]])){
  res_df <- temp[["imp_list"]][[cell_type]]$results
  res_df <- res_df[order(res_df$importance, decreasing = TRUE),]
  top_genes <- res_df$feature[1:top_n]

  genes_res_ls[[cell_type]] <- get_positive_contributors_for_class(
                                rf_model = rf_model_A,
                                class_name = cell_type,
                                data = rf_data, 
                                feature_to_fit = top_genes)
                              
}
save(genes_res_ls, file = "Rdata/10.hspc-combine.rf_model_A.iml-explanation-direction.genes_res_ls.rdata")

xlsx_file <- paste0("23.hspc.rf_model_A.iml.important-features-with-slopes.xlsx")
if ( file.exists(paste0("table/",xlsx_file)) ) {
  file.remove( paste0("table/",xlsx_file) )
}
wb <- my_xlsx_create(xlsx_file)

for (cell_type in names(genes_res_ls)){
  df <- genes_res_ls[[cell_type]]
  df$feature <- gsub("G_", "", df$feature)
  df$feature <- gsub("_", "-", df$feature)
  wb <- my_xlsx_write(wb, data = df, sheet = paste0(cell_type), rowNames = FALSE)
}
my_xlsx_save(wb, xlsx_file)
