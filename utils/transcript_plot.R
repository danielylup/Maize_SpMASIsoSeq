library(ggtranscript)
library(magrittr)
library(dplyr)
library(ggplot2)
library(rtracklayer)
library(stringr)
Zm_MAS_gtf <- rtracklayer::import("./GeneModel/Maize_LRiso.gtf")
VR03tleaf_IsoFraction_all <- read.delim("./VR03tleaf_IsoFraction_all.tsv")
colnames(VR03tleaf_IsoFraction_all) <- c("isoform_id", "IF", "gene_id", "PB_gene", "Read_count")

transcript_plot <- function(geneID, geneNAME){
  gtf <- as.data.frame(Zm_MAS_gtf) %>% 
    dplyr::filter(!is.na(gene_id), gene_id == geneID)
  exon <- gtf[gtf$type %in% "exon", c("seqnames", "start", "end", "width", "strand", "type", "gene_id", "gene_name", "transcript_id")]
  exon$transcript_id <- paste(exon$transcript_id, ":", sep = "")
  exon <- merge(exon, VR03tleaf_IsoFraction_all[, c("transcript_id", "IF")], by = "transcript_id")
  exon$gene_name <- paste(geneNAME, str_split_fixed(exon$transcript_id, "[.]", 3)[,3], sep = ".")
  exon %>%
    ggplot(aes(xstart = start, xend = end, y = gene_name)) +
        geom_range(aes(fill = IF)) +
        geom_intron(data = to_intron(exon, "transcript_id"), aes(strand = strand)) +
        ylab("Isoform") +
        theme_bw() +
        scale_fill_manual(values = c("red", "grey", "blue")) +
        theme(legend.key.size = unit(0.5, "cm"))
}
