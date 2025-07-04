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
library(cowplot) %>% suppressMessages()
library(sctransform) %>% suppressMessages()
library(stringr) %>% suppressMessages()
library(Matrix) %>% suppressMessages()
library(reshape2) %>% suppressMessages()
library(scales) %>% suppressMessages()
library(ggrepel) %>% suppressMessages()
library(ggh4x)  %>% suppressMessages()

### set path
setwd("/home/mayao/LABdata/MY/hsc1_hsc2/mpp")

### set Seurat
options(Seurat.object.assay.version = "v5")


### -----------------------
### 1. Calculate attributes
### -----------------------
##-----------------------##
load('Rdata/temp02.mpp.seurat.after_qc_norm.lognorm.rdata')
## reset seurat obj name as 'obj' 
obj <- mpp

##-----------------------##
## Gene attributes: mean and variance
## Cell attributes: total umi counts per cell, detected genes per cell
data <- obj[['RNA']]$scale.data
obj[['RNA']] <- JoinLayers(obj[['RNA']])

counts <- obj[['RNA']]$counts[rownames(data),colnames(data)]
meta <- obj@meta.data[colnames(counts),]

## Calculate gene attributes (raw umi counts)
ga <- calc_ga(counts)
message("#### ga ####")
print(head(ga))

## Calculate cell attributes
ca <- data.frame(row.names=rownames(meta), 
                    umi=meta$nCount_RNA, 
                    dg=meta$nFeature_RNA, 
                    mito=meta$pct_counts_Mito)
message("#### ca ####")
print(head(ca))

## Calculate gene attributes (normalized)
ga_norm <- calc_ga(data)
message("#### ga_norm ####")
print(head(ga_norm))


### --------------------
### 2. Metrics
### --------------------
# Things need to check:
# Metrics for cell quality 
# 1. Distribution of counts per cell
# 2. Distribution of detected genes per cell
# 3. Relationship between detected genes and counts
# 4. Relationship between mito.percent and counts

# Metrics for genes
# 5. Distribution of gene expression mean

# Metrics for gene attributes and counts relationship (needed for normalization)
# 6. Relationship between gene expression and counts (under a certain counts value [i.e., a cell], plot expression value of each gene)
# 7. Relationship between gene variance and counts (under a certain counts value [i.e., a cell], plot variance value of each gene)
# 8. Correlation of gene expression and cell umi counts (shown in gene groups)
# 9. Relationship between gene variance and gene mean (MA plot)

# Additional for SCT:
# 10. Distribution of corrected counts: nCount_SCT
# 11. Distribution of corrected features: nFeature_SCT

##-----------------------##
meta <- obj@meta.data %>% as.data.frame()

##-----------------------##
group_var <- 'Group'
is.umi <- TRUE # Refers to raw umi counts

##-----------------------##
## 1. Distribution of counts per cell
p1<- meta %>% 
        ggplot(aes(x=nCount_RNA, fill=!!ensym(group_var), color=!!ensym(group_var))) + 
        geom_density(alpha = 0.2) + 
        theme_classic() +
        scale_x_log10()
p1

## 2. Distribution of detected genes per cell
p2<- meta %>% 
        ggplot(aes(x=nFeature_RNA, fill=!!ensym(group_var), color=!!ensym(group_var))) + 
        geom_density(alpha = 0.2) + 
        theme_classic() +
        scale_x_log10()
p2

## 3. Relationship between detected genes and counts
p3 <- FeatureScatter(obj, 
                    feature1 = "nCount_RNA", 
                    feature2 = "nFeature_RNA",
                    group.by=group_var)
p3

# 4. Relationship between mito.percent and counts
p4 <- FeatureScatter(obj, 
                    feature1 = "nCount_RNA",
                    feature2 = "pct_counts_Mito",
                    group.by=group_var)
p4

##-----------------------##
cell_quality_pl <- c(list(p1),list(p2),list(p3),list(p4))
save_plotlist(cell_quality_pl, ncol=1, outfile=paste0("plot/QC.filtered_norm.1.cell_quality.pdf"), width=6, height=15)


##-----------------------##
# Put genes and cells into groups
ga <- g_in_groups(ga,4)[['ga']]
gene.grp.breaks <- g_in_groups(ga,4)[['gene.grp.breaks']]
gene.grp <- g_in_groups(ga,4)[['gene.grp']]

ca <- c_in_groups(ca,6)[['ca']]
cell.grp.breaks <- c_in_groups(ca,6)[['cell.grp.breaks']]
cell.grp <- c_in_groups(ca,6)[['cell.grp']]

##-----------------------##
# 5. Distribution of gene expression mean
p5 <- ggplot(ga, aes(x=log_gmean)) + geom_histogram(aes(fill=group), breaks=gene.grp.breaks) +
        scale_fill_manual(values = rev(brewer_pal(palette='YlOrRd')(length(levels(ga$group))+1)),
                            labels = table(ga$group),
                            name = 'Group size') +
        xlab('Gene mean [log10]') + ylab('Number of genes') +
        ggtitle('Histogram of mean gene UMI counts')
ga_pl <- c(list(p5))

##-----------------------##
save_plotlist(ga_pl, ncol=1, outfile=paste0("plot/QC.filtered_norm.2.gene_mean_distribution.pdf"), width=6, height=4)


##-----------------------##
trends <- get.trends(data = as.matrix(data), 
                        cell.umi = ca$umi, 
                        gene.grp = ga$group, 
                        cell.grp = ca$group, 
                        do.scale = FALSE)
trends$fitdf$y <- log10(pmax(0.0001, trends$fitdf$y))
# 6. Relationship between gene expression and counts
p6 <- plot.fit(trends$fitdf, is.umi = is.umi) +
        xlab('Total cell UMI [log10]') + 
        ylab('Gene UMI count [log10]') +
        ggtitle('Gene UMI count as function of cell UMI count\nMean and interquartile range per group')

# 7. Relationship between gene variance and counts
p7 <- plot.var(trends$vardf)

# 8. Correlation of gene expression and cell umi counts (shown in gene groups)
p8 <- plot.cor(trends$cordf)

# 9. Relationship between gene variance and gene mean (MA plot)
p9 <- ggplot(ga_norm, aes(mean, variance)) + geom_point(shape=16) + geom_abline(slope=1, color='deeppink') +
    xlab('Gene mean') + ylab('Gene variance') +
    ggtitle('Gene variance as function of gene mean')

##-----------------------##
ga_counts_pl <- c(list(p6),list(p7), list(p8), list(p9))    
save_plotlist(ga_counts_pl, ncol=1, outfile=paste0("plot/QC.filtered_norm.3.ga_counts_relation.pdf"), width=6, height=15)



##-----------------------##
# 10. Distribution of corrected counts: nCount_SCT
p10 <- meta %>% 
        ggplot(aes(x=nCount_SCT, fill=!!ensym(group_var), color=!!ensym(group_var))) + 
        geom_density(alpha = 0.2) + 
        theme_classic() +
        scale_x_log10()
p10

# 11. Distribution of corrected features: nFeature_SCT
p11 <- meta %>% 
        ggplot(aes(x=nFeature_SCT, fill=!!ensym(group_var), color=!!ensym(group_var))) + 
        geom_density(alpha = 0.2) + 
        theme_classic() +
        scale_x_log10()
p11

## 12. Relationship between nFeature_SCT and nCount_SCT
p12 <- FeatureScatter(obj, 
                    feature1 = "nCount_SCT", 
                    feature2 = "nFeature_SCT",
                    group.by=group_var)
p12


# 13. Relationship between mito.percent and nCount_SCT
p13 <- FeatureScatter(obj, 
                    feature1 = "nCount_SCT",
                    feature2 = "pct_counts_Mito",
                    group.by=group_var)
p13

sct_pl <- c(list(p10),list(p11),list(p12),list(p13))
save_plotlist(sct_pl, ncol=1, outfile=paste0("plot/QC.filtered_norm.4.sct_cell_quality.pdf"), width=6, height=15)



message("#### R04 done ####")

