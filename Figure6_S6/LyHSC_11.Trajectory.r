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

### set plot options
### Define color palettes and plot themes
colLibrary <- colorRampPalette(brewer.pal(n = 5, name = "Spectral"))(5)
colGEX <- c("grey85", brewer.pal(7, "Reds"))      # color for gene expression
colCcy <- c("black", "blue", "darkorange")        # color for cellcycle phase
plotTheme <- theme_classic(base_size = 14)
## change colGroup
colGroup <- c("HSC1"="black", "HSC2"="red", 
              "N1"="#FF80C0","N2"="#188386", "Ly-I"="#808080"
              )
show_col(colGroup)

## set other options after being decided
colCls <- c("#4189C9", "#D791BF","#FFCF3F")
show_col(colCls)
colState <- c("3" = "#CE95C2", "2" = "#4F67B0", "7" = "#07A6ED", "1" = "#74D109",
              "5" = "#D1272A", "6" = "#FBAA19")
show_col(colState)

## set python options and load modules
python_path <- '/home/mayao/.conda/envs/seurat/bin/'
use_python(python_path)
py_run_string("import numpy as np")
py_run_string("import matplotlib.pyplot as pl")

### set path
setwd("/home/mayao/LABdata/MY/hsc1_hsc2/mpp/")

### set Seurat
options(Seurat.object.assay.version = "v5")



## ------------------------------------------------------------------------------------
load("Rdata/RCombineHSC1.02.hsc-n1n2n3.combine-obj.rdata")
obj <- combine %>% subset(subset = Group != "HSC2")
obj$Group <- factor(obj$Group, levels = c("HSC1", "N1", "N2", "Ly-I" ))
obj@reductions %>% names()
table(obj$Group)

## ------------------------------------------------------------------------------------
## 1. Diffusion map
# obs <- c("Group", "seurat_clusters")
# group.by <- "Group"

# return_list <- setup_anndata(obj, obs=obs, pca_name='pca', umap_name='umap.pca', n_pcs=6, n_neighbors=30)
# sc <- return_list[['sc']]
# adata <- return_list[['adata']]

# sc$tl$diffmap(adata)
# oupDR <- py_to_r(adata$obsm['X_diffmap'])
# rownames(oupDR) <- colnames(obj)
# colnames(oupDR) <- paste0("DC_", 0:(15-1))
# oupDR <- oupDR[, paste0("DC_", seq(3))]
# obj[["diffmap"]] = CreateDimReducObject(embeddings = oupDR, key = "DC_",
#                                         assay = DefaultAssay(obj))

# p3 <- dimplot(obj, reduction = "diffmap", group.by='Group', pt.size = 1.5, 
#             label = FALSE,
#             dims = c(1,2))+
#             force_panelsizes(rows=unit(6,"cm"), cols=unit(5,"cm"))
# p4 <- dimplot(obj, reduction = "diffmap", group.by='Group', pt.size = 1.5, 
#             label = FALSE,
#             dims = c(1,3))+
#             force_panelsizes(rows=unit(6,"cm"), cols=unit(5,"cm"))
# p3
# p4
# save_plotlist(plotlist=c(list(p3),list(p4)), ncol=2, outfile="plot/91.RCombineHSC1.hsc1hsc2-n1n2n3.diffusionmap.pdf", width=15, height=10) 


# library(plotly)
# dc <- Embeddings(obj, reduction = "diffmap") %>% as.data.frame()
# dc$Group <- obj$Group

# p5 <- plot_ly(x=as.numeric(dc$DC_1), y=dc$DC_2, z=dc$DC_3, type="scatter3d", mode="markers", color=dc$Group, colors = colGroup)
# htmlwidgets::saveWidget(p5, file = "plot/93.RCombineHSC1.hsc1hsc2-n1n2n3.3D-diffusionmap.html")


## Diffusion map
obs <- c("Group", "Phase")
group.by <- "Group"
pca_name <- "harmony"


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
htmlwidgets::saveWidget(p, file = paste0("plot/94.RCombineHSC1.hsc1hsc2-n1n2n3.3D-diffusionmap.",pca_name,".html"))




## -------------------------------------------------------------------------------------
## 3. Diffusion Map and PAGA
message('#### Diffusion Map and PAGA ####')
obs <- c("Group", "Phase")
group.by <- 'Group'
pca_name <- "harmony"
umap_name=paste0('umap.', pca_name)
nPC <- 6
colors <- colGroup

return_list <- setup_anndata(obj, obs=obs, pca_name=pca_name, umap_name=umap_name, n_pcs=nPC, n_neighbors=30)
sc <- return_list[['sc']]
adata <- return_list[['adata']]

sc$tl$diffmap(adata)
adata$obs

plotname1 <- paste0('plot/temp.hsc1hsc2-n1n2n3.PAGA.PAGAgraph.pdf')
plotname2 <- paste0('plot/temp.hsc1hsc2-n1n2n3.PAGA.Group.pdf')
plotname3 <- paste0('plot/temp.hsc1hsc2-n1n2n3.PAGA.dpt.pdf')
# Run paga
sc$tl$paga(adata, groups = group.by)

# plot paga graph G*
sc$pl$paga(adata, color= group.by, cmap = colors)
py_run_string("fig1 = pl.gcf()")
py_run_string("fig1.savefig(r.plotname1)")

# initializing using paga
sc$tl$draw_graph(adata, init_pos='paga')

# draw paga in single-cell resolution 
sc$pl$draw_graph(adata, color=group.by, palette = colors, return_fig=TRUE)
py_run_string("fig2 = pl.gcf()")
py_run_string("fig2.savefig(r.plotname2)")

# run dpt
py_run_string(paste0("r.adata.uns['iroot'] = np.flatnonzero(r.adata.obs['Group']  == 'HSC1')[0]"))
sc$tl$dpt(adata)

# draw dpt
sc$pl$draw_graph(adata, color=c(group.by, 'dpt_pseudotime'), legend_loc='on data', return_fig=TRUE)
py_run_string("fig3 = pl.gcf()")
py_run_string("fig3.savefig(r.plotname3)")

# add to seurat meta
obj$paga_dpt <- py_to_r(adata$obs$dpt_pseudotime)
obj$paga_rank <- rank(obj$paga_dpt)

message("**** Plot boxplot ****")
p <- pseudotime_boxplot(obj, group.by='Group', rank.by='paga_rank', group.colors=colCls)
p




## -------------------------------------------------------------------------------------
## 4 Monocle
active_idents <- 'Group'
reduc <- 'umap.harmony'
rootCellType <- 'HSC1'
assay <- "RNA"
cols_div <- c("#EDD1CB", "#DAA4AC", "#BD7A98", "#935685", "#613969", "#2D1E3E") 
# cols_map <- colorRampPalette(cols_div)(100)
cols_map <- viridisLite::mako(100)
obs <- c("Group", "Phase")

## Set idents
Idents(obj) <- active_idents

## Get cds
cds <- seurat_to_monocle(obj, assay=assay, embeddings=reduc)

## Learn trajectory 
cds <- learn_graph(cds, use_partition = F)

# ## Order
# cds <- order_cells(cds, reduction_method='UMAP', root_cells = colnames(cds[, clusters(cds) == rootCellType]))

# ## Plot pseudotime trajectory
# # plotcol <- cols_map[cut(pseudotime(cds), breaks=100)]
# my_plot_cells(cds, 
#                 reduction_method = 'UMAP',
#                 color_cells_by = "cluster", 
#                 label_groups_by_cluster = F,
#                 label_branch_points = F, 
#                 label_roots = FALSE, 
#                 label_leaves = FALSE,
#                 cell_size = 1.5,
#                 alpha = 0.5,
#                 plotcol_div=colGroup,
#                 # plotcol_map=cols_map,
#                 trajectory_graph_color="black",
#                 trajectory_graph_segment_size=0.8,
#                 width = 6,
#                 height = 5)

## get principal graph
principal_coords <- cds@principal_graph_aux[["UMAP"]]$dp_mst %>% t() %>% as.data.frame()
head(principal_coords)
graph_edges <- principal_graph(cds)$UMAP %>% igraph::as_data_frame(what = "edges")
head(graph_edges)

## plot 3D UMAP
library(plotly)
umap_coords <- reducedDims(cds)$UMAP %>% as.data.frame()
umap_coords$Group <- obj$Group

p <- plot_ly(x=as.numeric(umap_coords$umapharmony_1), y=umap_coords$umapharmony_2, z=umap_coords$umapharmony_3, type="scatter3d", mode="markers", color=umap_coords$Group, colors = colGroup)

for (i in 1:nrow(graph_edges)) {
  from <- graph_edges$from[i]
  to <- graph_edges$to[i]

  x_line <- c(principal_coords[from, 1], principal_coords[to, 1])
  y_line <- c(principal_coords[from, 2], principal_coords[to, 2])
  z_line <- c(principal_coords[from, 3], principal_coords[to, 3])

  p <- add_trace(p,
    x = x_line, y = y_line, z = z_line,
    type = 'scatter3d', mode = 'lines',
    line = list(color = 'black', width = 3),
    showlegend = FALSE,
    inherit = FALSE
  )
}

htmlwidgets::saveWidget(p, file = paste0("plot/97.RCombineHSC1.hsc1hsc2-n1n2n3.3D-umapharmony-Monocle3.html"))



# plot_genes_in_pseudotime(cds[c("Gata2", "Runx1", "Itgam"),])
# ggsave("plot/84.RCombineHSC1.hsc1-n1n2n3.monocle3.pdf")
ggsave("plot/92.RCombineHSC1.hsc1hsc2-n1n2n3.monocle3.pdf")
