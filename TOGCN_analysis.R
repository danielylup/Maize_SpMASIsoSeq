#TO-GCN analysis requires normalized data (TPM) in time-series dataset
library(GenomicFeatures)
library(dplyr)
library(stringr)

# Import read count matrix from HTseq 
MaizeV5ext_HTseq <- read.delim("./LCM_data/HTSeq_output/Maize_V5ext_HTseq.txt", header=FALSE)

# Proces read count matrix and normalizes using TPM
rownames(MaizeV5ext_HTseq) <- str_split_fixed(MaizeV5ext$V1, "_PB", 2)[,1]
MaizeV5ext_HTseq <- MaizeV5ext_HTseq[,2:79]

# Retrieve exonic length of each gene
MaizeV5ext_GFF <- makeTxDbFromGFF("./gene_model/Maize_V5ext.gff3", format = "gff3")
exon_list_per_gene <- exonsBy(MaizeV5ext_GFF, by="gene")
exonic_gene_sizes <-sum(width(reduce(exon_list_per_gene)))
MaizeV5ext_exonic_length <- (as.data.frame(exonic_gene_sizes))
rownames(MaizeV5ext_exonic_length) <- make.names(str_split_fixed(rownames(MaizeV5ext_exonic_length), "_PB", 2)[,1], unique = TRUE)
# Compile the read count matrix data with the exonic length column
MaizeV5ext_TPM <- merge(MaizeV5ext_HTseq, MaizeV5ext_exonic_length, by = 0)

# Function to calculate TPM
tpm_mat <- function(counts,length) {
    x <- counts/length
    return(t(t(x)*1000000/sum(x)))
    }

# calculate the TPM of each gene across all LCM RNA-seq datasets
for(i in 2:((ncol(MaizeV5ext_TPM))-1)){
    MaizeV5ext_TPM[,i] <- tpm_mat(MaizeV5ext_TPM[,i], MaizeV5ext_TPM$exonic_gene_sizes)
    }
# Data frame processing
rownames(MaizeV5ext_TPM) <- MaizeV5extV6_GeneModel_TPM$Row.names
MaizeV5ext_TPM <- MaizeV5ext_TPM[, -1]

# Import the metadata for all the LCM RNA-seq dataset
SRR_LCM_data <- read.delim("./LCM_data/SRR_LCM_data.txt", header = TRUE)
# Insert required information from metadata into the TPM data frame
MaizeV5ext_TPM_annot <- as.data.frame(merge(t(MaizeV5ext_TPM), SRR_LCM_data[,1:2], by.x = 0, by.y = 'SRA'))
MaizeV5ext_TPM_annot$stages <- str_split_fixed(MaizeV5ext_TPM_annot$Library, "_", 2)[,1]
MaizeV5ext_TPM_aggr <- MaizeV5ext_TPM_annot %>%
    group_by(stages) %>%
    summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE))) %>%
    as.data.frame()
rownames(MaizeV5ext_TPM_aggr) <- MaizeV5ext_TPM_aggr$stages
MaizeV5ext_TPM_aggr <- MaizeV5ext_TPM_aggr[,-1]
MaizeV5ext_TPM_aggr <- as.data.frame(t(MaizeV5ext_TPM_aggr[c("mGM", "3C", "4BS+V", "5BS+V", "6BS+V", "PM", "3PM", "4PM", "5/6PM", "1M", "2M"), ]))

#Prepare TF expression list and rename TF gene
MaizeV5ext_TPM_TF <- MaizeV5ext_TPM_aggr[rownames(MaizeV5ext_TPM_aggr) %in% ZmTF_UniqList$ZmV5 ,]
rownames(MaizeV5ext_TPM_TF) <- paste(rownames(MaizeV5ext_TPM_TF), "TF", sep = "_")

MaizeV5ext_TPM_aggr <- MaizeV5ext_TPM_aggr[!rownames(MaizeV5ext_TPM_aggr) %in% ZmTF_UniqList$ZmV5 ,]
MaizeV5ext_TPM_aggr <- rbind(MaizeV5ext_TPM_aggr, MaizeV5ext_TPM_TF)

#exprt tsv file of both all gene and TF gene files
write.table(MaizeV5ext_TPM_aggr, "./TOGCN_analysis/MaizeV5ext_TPM_aggr.tsv", sep = '\t')
write.table(MaizeV5ext_TPM_TF, "./TOGCN_analysis/MaizeV5ext_TPM_TF.tsv", sep = '\t')

########################################################################################################################################################################################
TO-GCN run in linux

Cutoff 5 ./TOGCN_analysis/MaizeV5ext_TPM_TF.tsv 
TO-GCN 5 ./TOGCN_analysis/MaizeV5ext_TPM_TF.tsv ./TOGCN_analysis/seeds.txt 0.92   # Zm00001eb363810 TF used as the seed that should have high expression in early stage
GeneLevel 5 ./TOGCN_analysis/MaizeV5ext_TPM_TF.tsv ./TOGCN_analysis/MaizeV5ext_TPM_aggr.tsv ./TOGCN_analysis/Node_level.tsv 0.92

# Visualization of the TO-GCN result is done using Cytoscape, with "Node_relation.csv" as main input and "Node_level.tsv" to arrange them based on TO-GCN level

########################################################################################################################################################################################
# Generate heatmap of TF expression across TO-GCN level
library(dplyr)
library(pheatmap)

TFlevel_MaizeV5ext <- read.csv("./TOGCN_analysis/Gene_list_in_each_level.csv")
colnames(TFlevel_MaizeV5ext) <- c("Gene_ID", "is_TF", "TOGCN_level")

# Generate data frame that consist of TF's mean expression in each TO-GCN level and cell types
TF_list <- list()
for (i in 1:tail(unique(TFlevel_MaizeV5ext$TOGCN_level), 1)){
  TF_temp <- colMeans(MaizeV5ext_TPM_TF[rownames(MaizeV5ext_TPM_TF) %in% (subset(TFlevel_MaizeV5ext, is_TF == "1" & TOGCN_level == i))$Gene_ID,])
  TF_temp <- TF_temp[c("mGM", "3C", "4BS+V", "5BS+V", "6BS+V")]
  TF_list <- append(TF_list, list(TF_temp))
}
# Reorganize the data frame to use for plotting
TOGCN_TFmat <- TF_list %>%
  bind_rows() %>%
  scale() %>%
  as.data.frame()
# Append "L" in front to significant the word "Level"
rownames(TOGCN_TFmat) <- paste("L", 1:tail(unique(TFlevel_MaizeV5ext$TOGCN_level), 1), sep = "")
# Generating the heatmap plot
pheatmap(TOGCN_TFmat, cluster_cols = F, cluster_rows = F, angle_col = 0, scale = 'row', display_numbers = TRUE)


#######################################################################################################################################################################
library(ggsankey)
# The exact same TO-GCN analysis was conducted using the gene expression matrix from Maize_V5can gene model as reference

# Construct data frame for sankey plot
TFlist_level_compare <- merge(subset(TFlevel_MaizeV5can, is_TF == "1"), subset(TFlevel_MaizeV5ext, is_TF == "1"), by = "Gene_ID")
TFlist_level_sankey <- TFlist_level_compare %>%
  ggsankey::make_long("TOGCN_level.x", "TOGCN_level.y")

# Sankey plot compares the TF's TO-GCN level assignment between Maize_V5can and Maize_V5ext gene models
ggplot(TFlist_level_sankey, aes(x=x, next_x = next_x, node = node, next_node =next_node, fill = factor(node), label = node)) +
    geom_sankey(flow.alpha = 0.5) +
    geom_sankey_label(size = 3.5, color = 1, fill = "white") +
    scale_fill_viridis_d(alpha = 0.8) +
    theme_sankey(base_size =  16) +
    xlab("Gene model") +
    theme(legend.position = "none") +
    coord_flip()
