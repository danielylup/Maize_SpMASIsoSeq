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



##### Pseudotime Trajectory Analysis #####
library(monocle3)
library(SeuratWrappers)


cds <- as.cell_data_set(VR03tleaf_WGCNA)
# Generate trajectory on UMAP 
cds <- cluster_cells(cds, reduction_method='UMAP')
cds <- learn_graph(cds, use_partition = FALSE, learn_graph_control = list(prune_graph=TRUE, ncenter=500))

VR03tleaf_seurat_WGCNA[["UMAP"]] <- VR03tleaf_seurat_WGCNA[["umap_iso"]]
# convert the seurat object to CDS
cds <- as.cell_data_set(VR03tleaf_seurat_WGCNA)
# run the monocle clustering
cds <- cluster_cells(cds, reduction_method='UMAP')
# learn graph for pseudotime
cds <- learn_graph(cds)

# plot the trajectory graph on UMAP plot
p1 <- plot_cells(
  cds = cds,
  color_cells_by = "embryonic.vein",
  show_trajectory_graph = TRUE,
  label_principal_points = TRUE,
  cell_size = 1.5,
  group_label_size = 5
) 

# Determined principal node based on the biology (SAM) and CytoTRACE analysis
principal_node <- 'Y_5'
cds <- order_cells(cds,root_pr_nodes = principal_node)

VR03tleaf_WGCNA$pseudotime <- pseudotime(cds)
VR03tleaf_WGCNA$vein_pseudotime <- ifelse(VR03tleaf_WGCNA$embryonic.vein %in% c("SAM_P1", "P2_P3", "embry_v1", "embry_v2"), VR03tleaf_WGCNA$pseudotime, NA)
VR03tleaf_WGCNA$leaf_pseudotime <- ifelse(VR03tleaf_WGCNA$embryonic.vein %in% c("SAM_P1", "P2_P3", "P4_P5", "shoot_col"), VR03tleaf_WGCNA$pseudotime, NA)

# intergrate UMAP information into meta.data for plotting
VR03tleaf_WGCNA$UMAP1 <- VR03tleaf_WGCNA@reductions$umap_iso@cell.embeddings[,1]
VR03tleaf_WGCNA$UMAP2 <- VR03tleaf_WGCNA@reductions$umap_iso@cell.embeddings[,2]

pseu_leaf <- VR03tleaf_seurat_WGCNA@meta.data %>%
    ggplot(aes(x=UMAP1, y=UMAP2, color=leaf_pseudotime)) +
    ggrastr::rasterise(geom_point(size=3), dpi=500, scale=0.75) +
    coord_equal() +
    scale_color_gradientn(colors=viridis(256), na.value='grey') +
    umap_theme()

pseu_vein <- VR03tleaf_seurat_WGCNA@meta.data %>%
    ggplot(aes(x=UMAP1, y=UMAP2, color=vein_pseudotime)) +
    ggrastr::rasterise(geom_point(size=3), dpi=500, scale=0.75) +
    coord_equal() +
    scale_color_gradientn(colors=plasma(256), na.value='grey') +
    umap_theme()

#plot module dynamics in pseudotime trajectory
PlotModuleTrajectory(
    VR03tleaf_seurat_WGCNA,
    pseudotime_col = 'leaf_pseudotime'
)

PlotModuleTrajectory(
    VR03tleaf_seurat_WGCNA,
    pseudotime_col = 'vein_pseudotime'
)



##### Functional Enrichment Analysis #####
# Functional enrichment using the top 25 hub isoforms for each co-expression modules
library(gprofiler2)
library(rtracklayer)

# Create a list of PacBio isoform ID to indentify the isoform
LRiso_gtf <- readGFF("./gene_model/LRiso_GeneModel.gtf")
PB_newname_list <- data.frame()[FALSE,"Name"]
PB_newname_list$Name <- LRiso_gtf[grep("transcript", LRiso_gtf$type), "transcript_id"]
PB_newname_list$Name <- gsub("XLOC_", "XLOC", PB_newname_list$Name)
PB_newname_list$Name <- paste(PB_newname_list$Name, ":", sep = "")
PB_newname_list$Isoform <- str_split_fixed(PB_newname_list$Name, "_", 2)[,2]
PB_newname_list$Gene <- str_split_fixed(PB_newname_list$Name, "_", 2)[,1]
PB_newname_list <- as.data.frame(PB_newname_list)

# Extract the hub isofomr for each co-expression modules
VR03tleaf_hubgene <- GetHubGenes(VR03tleaf_WGCNA, n_hubs = 25)
colnames(VR03tleaf_hubgene) <- c("Isoform", "module", "kME")
VR03tleaf_hubgene <- (merge(VR03tleaf_hubgene, PB_newname_list, by = "Isoform"))

# GO term enrichment using gprofiler2
for(i in 1:length(mod)) {
  filenames <-paste0("./hdWGCNA_analysis/Module_gprofiler_tleaf_25hubg/module_", mod[i], ".png")
  VR03tleaf_gost <- gost((subset(VR03tleaf_hubgene, module == mod[i]))$Gene, organism = "zmays", sources = c("GO:BP", "KEGG"), highlight = TRUE)
  if(is.null(VR03tleaf_gost)){
    print(paste0("There is no enriched GO term in module ", mod[i], "! Skip output table!"))
  } else {
  publish_gosttable(VR03tleaf_gost, highlight_terms = subset(VR03tleaf_gost$result, highlighted == "TRUE"), 
                    show_columns = c("term_id", "term_name", "term_size", "intersection_size", "p_value"), 
                    filename = filenames)
                    }
}

# Due to the output is only available in png format, I manually input all the enriched GO terms into dataset format
VR03tleaf_25hubg_GOBP <- read.csv("./hdWGCNA_analysis/VR03tleaf_25hubg_GOBP.csv")

# Plot the GO term enrichment bubble plots
ggplot(VR03tleaf_25hubg_GOBP, aes(x=module, y=reorder(term_name, sorting), colour=`p.value`, size=enrichment_ratio)) +
    geom_point() +
    scale_color_gradient(high = "mediumblue", low = "red2") +
    expand_limits(x=0) +
    labs(x="Enrichment ratio", y="GO term (biological process)", colour="p-value", size=paste('Enrichment\nratio')) +
    theme_bw() +
    theme(legend.position = c(0.8, 0.25), legend.background = element_rect(colour = "black"), axis.text = element_text(colour = "black", size = 10, family = "arial"))
