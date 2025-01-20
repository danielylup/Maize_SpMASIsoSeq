library(Seurat)
library(ggplot2)
library(ggrepel)
library(WGCNA)
library(hdWGCNA)
library(SeuratWrappers)
library(patchwork)
library(dplyr)
library(gganimate)
library(harmony)

DefaultAssay(VR03tleaf_seurat) <- 'iso'
VR03tleaf_WGCNA <- SetupForWGCNA(VR03tleaf_seurat,  wgcna_name = "tleaf")
# construct metacells:
VR03tleaf_WGCNA <- MetacellsByGroups(VR03tleaf_WGCNA, 
                                     group.by = 'embryonic.vein', 
                                     ident.group = "embryonic.vein", 
                                     k = 50, 
                                     max_shared=25, 
                                     min_cells=75, 
                                     reduction = 'harmony_iso', 
                                     assay = 'iso', 
                                     mode = 'sum')
VR03tleaf_WGCNA <- NormalizeMetacells(VR03tleaf_WGCNA)

groups <- unique(GetMetacellObject(VR03tleaf_WGCNA)$embryonic.vein)
# set up gene expression matrix
VR03tleaf_WGCNA <- SetDatExpr(VR03tleaf_WGCNA, 
                              group.by='embryonic.vein', 
                              group_name = groups, 
                              use_metacells=TRUE, 
                              slot = 'data', 
                              assay = 'iso')
# test soft powers
VR03tleaf_WGCNA <- TestSoftPowers(VR03tleaf_WGCNA)
plot_list <- PlotSoftPowers(VR03tleaf_WGCNA)
wrap_plots(plot_list, ncol=2)

# construct wgcna network (auto-detect soft_power):
VR03tleaf_WGCNA <- ConstructNetwork(VR03tleaf_WGCNA, 
                                    setDatExpr=FALSE, 
                                    detectCutHeight=0.980, 
                                    mergeCutHeight=0.2, 
                                    TOM_name = 'tleaf', 
                                    minModuleSize=100, 
                                    overwrite_tom = TRUE)

# Plot dendrogram
PlotDendrogram(VR03tleaf_WGCNA, main='hdWGCNA Dendrogram')

# compute module connectivity:
VR03tleaf_WGCNA <- ModuleEigengenes(VR03tleaf_WGCNA, verbose = TRUE)
modules_orig <- GetModules(VR03tleaf_WGCNA)
VR03tleaf_WGCNA <- ModuleConnectivity(VR03tleaf_WGCNA)

# get the MEisos from the seurat object and add it to the metadata
MEiso <- GetMEs(VR03tleaf_WGCNA)
meta <- VR03tleaf_WGCNA@meta.data
VR03tleaf_WGCNA@meta.data <- cbind(meta, MEiso)

# get a list of features to plot
modules <- GetModules(VR03tleaf_WGCNA)
mods <- levels(modules$module)
mods <- mods[mods!='grey']

# make DotPlot
DotPlot(VR03tleaf_WGCNA,
        group.by='embryonic.vein',
        features = rev(mods)
) + RotatedAxis() +
  scale_color_gradient2(high='red', mid='grey95', low='blue') + xlab('') + ylab('') +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1)
  )

# Plot UMAP connectivity plot
VR03tleaf_WGCNA <- RunModuleUMAP(VR03tleaf_WGCNA, n_hubs =25, n_neighbors=15, min_dist=0.2, spread=1)
ModuleUMAPPlot(VR03tleaf_WGCNA, edge.alpha=0.5, sample_edges=TRUE, keep_grey_edges=FALSE, edge_prop=0.075)

#plot co-expression network for each module
ModuleNetworkPlot(VR03tleaf_WGCNA)
