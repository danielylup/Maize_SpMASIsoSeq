###### EviAnn - evidence-based eukaryotic genome annotation software ######
# EviAnn run using uncollapsed FLNC transcripts sequences fastq file as input, and whole genome reference sequence and protein reference sequences retrieved from NCBI as reference input (-g & -r, respectively)
eviann.sh -t 56 \
          -g ./Zea_mays.Zm-B73-REFERENCE-NAM-5.0.genomics.fna \
          -e ./VR03.Skera.s_reads.segmented.fltnc.corrected.sorted.dedup.ensemble.mapped.fastq \
          -r ./Maize_V5ext/EviAnn_analysis/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.protein.fa

# Identify gene models that are also annotated from EviAnn analysis using gffcompare
gffcompare -r ./Maize_V5ext/EviAnn_analysis/ZmEviAnn_output.gtf \
           -s ./Zea_mays.Zm-B73-REFERENCE-NAM-5.0.genomics.fna \
           -o ./Maize_V5ext/EviAnn_analysis/gffcompare_Novel_EviAnn \
           ./gene_model/LRiso_GeneModel.gtf

# Calculate novel gene model with support from EviAnn annotation
grep "XLOC_" gffcompare_Novel_EviAnn.Zea_mays-B73-REFERENCE_V5ext_sorted.gtf.tmap | 
  awk '{if($1 != "-" && $3 != "u" && $3 != "i" && $3 != "y" && $3 != "s" && $3 != "x") print}' | 
  grep -v "Zm00001eb" | 
  cut -f4 | sort | uniq | wc -l
