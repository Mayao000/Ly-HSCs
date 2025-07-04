# '''
# This is to analyse pdc1 and pdc2.
# '''
rm(list=ls())
set.seed(1234)

print('---------------------------------')
print(Sys.Date())
print('---------------------------------')


source("/home/mayao/script/Function.r")

### ------------------------------------------------
###                   Running Zone
### ------------------------------------------------

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
# library(harmony) %>% suppressMessages()
library(scales) %>% suppressMessages()
### set path
setwd("/home/mayao/LABdata/MY/pdc")
### set Seurat
options(Seurat.object.assay.version = "v5")



### -----------------------
### 1. Calculate attributes
### -----------------------
##-----------------------##
load("Rdata/temp01.pdc.seuratObj.rdata" )
## reset seurat obj name as 'obj' 
obj <- pdc

##-----------------------##
## Gene attributes: mean and variance
## Cell attributes: total umi counts per cell, detected genes per cell
ga_ls <- list()
ca_ls <- list()
for(layer_n in Layers(obj)){
    counts <- obj$RNA[layer_n]
    meta <- obj@meta.data[colnames(counts),]

    ga <- calc_ga(counts)
    print(head(ga))
    ga_ls <- c(ga_ls, list(ga))

    ca <- data.frame(row.names=rownames(meta), 
                     umi=meta$nCount_RNA, 
                     dg=meta$nFeature_RNA, 
                     mito=meta$pct_counts_Mito)
    print(head(ca))
    ca_ls <- c(ca_ls, list(ca))
}
names(ga_ls) <- gsub('counts.','',Layers(obj))
names(ca_ls) <- gsub('counts.','',Layers(obj))

### --------------------
### 2. Before filtering
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

##-----------------------##
meta <- obj@meta.data %>% as.data.frame()

##-----------------------##
group_var <- "Name"
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
p1+p2

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
save_plotlist(cell_quality_pl, ncol=1, outfile=paste0("plot/QC.before_filter_raw.1.cell_quality.pdf"), width=6, height=15)


##-----------------------##
for (i in 1:length(ga_ls)){
    ga <- ga_ls[[i]]
    # remove genes whose logmean is NA
    ga <- ga[ga$log_gmean!=-Inf,]

    ca <- ca_ls[[i]]
    layer_n <- names(ga_ls)[i]
    counts <- obj$RNA[layer_n]

    ##-----------------------##
    # Put genes and cells into groups
    ga <- g_in_groups(ga,6)[['ga']]
    gene.grp.breaks <- g_in_groups(ga,6)[['gene.grp.breaks']]
    gene.grp <- g_in_groups(ga,6)[['gene.grp']]

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
    save_plotlist(ga_pl, ncol=1, outfile=paste0("plot/QC.before_filter_raw.2.gene_mean_distribution.",names(ga_ls)[i],".pdf"), width=6, height=4)


    ##-----------------------##
    trends <- get.trends(data = as.matrix(counts), 
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
    p9 <- ggplot(ga, aes(log_gmean, log10(variance))) + geom_point(shape=16) + geom_abline(slope=1, color='deeppink') +
    xlab('Gene mean [log10]') + ylab('Gene variance [log10]') +
    ggtitle('Gene variance as function of gene mean')

    ##-----------------------##
    ga_counts_pl <- c(list(p6),list(p7), list(p8), list(p9))    
    save_plotlist(ga_counts_pl, ncol=1, outfile=paste0("plot/QC.before_filter_raw.3.ga_counts_relation.",names(ga_ls)[i],".pdf"), width=6, height=15)
}
message("#### R02 done ####")

