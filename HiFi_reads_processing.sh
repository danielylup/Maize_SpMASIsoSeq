## HiFi reads processing and isoform structure classification using IsoSeq3 and Pigeon pipeline

# HiFi read segmentation: 
skera split -j 56 \
            --log-level INFO \
	    QD003_m64531e_230118_204124.hifi_reads.bam \
            mas16_primers.fasta \
	    QD003.Skera.s_reads.segmented.bam
skera split -j 56 \
            --log-level INFO \
	    QD010_m64531e_230415_193517.hifi_reads.bam \
            mas16_primers.fasta \
	    QD010.Skera.s_reads.segmented.bam
skera split -j 56 \
            --log-level INFO \
	    QD011_m64531e_230421_082528.hifi_reads.bam \
            mas16_primers.fasta \
	    QD011.Skera.s_reads.segmented.bam
# Skera version 0.1.0  
# Thread number setting available (-j [default 0])
# Output: *.Skera.s_reads.segmented.bam
# Summary output: *.Skera.s_reads.segmented_2.summary.csv
# Approximate run time: 4m 23s (56 threads)

#Primer removal:
lima -j 56 \
     --log-level INFO \
     QD003.Skera.s_reads.segmented.bam \
     10x_3kit_primers.fasta \
     QD003.Skera.s_reads.segmented.primer_rm.bam \
     --isoseq
lima -j 56 \
     --log-level INFO \
     QD010.Skera.s_reads.segmented.bam \
     10x_3kit_primers.fasta \
     QD010.Skera.s_reads.segmented.primer_rm.bam \
     --isoseq
lima -j 56 \
     --log-level INFO \
     QD011.Skera.s_reads.segmented.bam \
     10x_3kit_primers.fasta \
     QD011.Skera.s_reads.segmented.primer_rm.bam \
     --isoseq
# lima version 2.6.99 
# Thread number setting available (-j [default 0])
# Output: *.Skera.s_reads.segmented.primer_rm.bam		
# Summary output: *.Skera.s_reads.segmented.primer_rm.lima.summary
# Approximate run time: 8m 22 (56 threads)

# Tag clipped:
isoseq3 tag -j 56 \
            --log-level INFO \
	    QD003.Skera.s_reads.segmented.primer_rm.bam \
            QD003.Skera.s_reads.segmented.flt.bam \
	    --design T-12U-16B
isoseq3 tag -j 56 \
            --log-level INFO \
	    QD010.Skera.s_reads.segmented.primer_rm.bam \
            QD010.Skera.s_reads.segmented.flt.bam \
	    --design T-12U-16B
isoseq3 tag -j 56 \
            --log-level INFO \
	    QD011.Skera.s_reads.segmented.primer_rm.bam \
            QD011.Skera.s_reads.segmented.flt.bam \
	    --design T-12U-16B
#isoseq3 version 3.8.1 
#thread number setting available (-j [default 0])
#output: *.Skera.s_reads.segmented.flt.bam
#Summary output: None
#Approximate run time: 6m 38s (56 threads)

#Refine:
isoseq3 refine -j 56 \
               --log-level INFO \
	       --require-polya \
	       QD003.Skera.s_reads.segmented.flt.bam \
	       10x_3kit_primers.fasta \
	       QD003.Skera.s_reads.segmented.fltnc.bam \
isoseq3 refine -j 56 \
               --log-level INFO \
	       --require-polya \
	       QD010.Skera.s_reads.segmented.flt.bam \
	       10x_3kit_primers.fasta \
	       QD010.Skera.s_reads.segmented.fltnc.bam \
isoseq3 refine -j 56 \
               --log-level INFO \
	       --require-polya \
	       QD011.Skera.s_reads.segmented.flt.bam \
	       10x_3kit_primers.fasta \
	       QD011.Skera.s_reads.segmented.fltnc.bam \
# Thread number setting available (-j [default 0])
# Output: *.Skera.s_reads.segmented.fltnc.bam
# Summary output: *.Skera.s_reads.segmented.fltnc.report.csv
# Approximate run time: 3m 34s (56 threads)

# Merge all BAM files 
ls QD003.Skera.s_reads.segmented.fltnc.bam QD010.Skera.s_reads.segmented.fltnc.bam QD011.Skera.s_reads.segmented.fltnc.bam > VR03.Skera.s_reads.segmented.fltnc.fofn

# Cell barcode correction:
isoseq3 correct -j 56 \
                --log-level INFO \
		--barcodes visium_v1.RC.txt \
                VR03.Skera.s_reads.segmented.fltnc.bam \
		VR03.Skera.s_reads.segmented.fltnc.corrected.bam
# Thread number setting available (-j [default 0])
# Output: VR03.Skera.s_reads.segmented.fltnc.corrected.bam
# Summary output: VR03.Skera.s_reads.segmented.fltnc.corrected.report.json
# Approximate run time: 335s (56 threads)

# Deduplication:	   
samtools sort -@ 56 \
              -t CB VR03.Skera.s_reads.segmented.fltnc.corrected.bam \
	      -o VR03.Skera.s_reads.segmented.fltnc.corrected.sorted.bam
isoseq3 groupdedup -j 56 \
                   --log-level INFO \ 
		   VR03.Skera.s_reads.segmented.fltnc.corrected.sorted.bam \
                   VR03.Skera.s_reads.segmented.fltnc.corrected.sorted.dedup.bam
# Thread number setting available (-@ & -j [default 0])
# Output: VR03.Skera.s_reads.segmented.fltnc.corrected.sorted.dedup.bam
# Summary output: VR03.Skera.s_reads.segmented.fltnc.corrected.report.json
# Approximate run time: 22m 27s (56 threads; not including samtools)

# Map reads to a reference genome:    
pbmm2 align -j 56 \
            --log-level INFO \
	    --preset ISOSEQ \
            --sort VR03.Skera.s_reads.segmented.fltnc.corrected.sorted.dedup.bam \
	    Zea_mays.Zm-B73-REFERENCE-NAM-5.0.dna.toplevel.fa \
            VR03.Skera.s_reads.segmented.fltnc.corrected.sorted.dedup.ensemble.mapped.bam 
# pbmm2 version 1.9.0
# Output: VR03.Skera.s_reads.segmented.fltnc.corrected.sorted.dedup.ensemble.mapped.bam
# Summary output:
# Approximate run time:

# Collapse into unique isoforms:	   
isoseq3 collapse -j 56 \
		 --log-level INFO \
                 VR03.Skera.s_reads.segmented.fltnc.corrected.sorted.dedup.ensemble.mapped.bam \
		 VR03.Skera.s_reads.segmented.fltnc.corrected.sorted.dedup.ensemble.mapped.collapsed.gff
# Output: VR03.Skera.s_reads.segmented.fltnc.corrected.sorted.dedup.ensemble.mapped.collapsed.gff
# Summary output: 
# Approximate run time: 3m 27s

# Sort input transcript GFF: 	   
pigeon sort -j 56 \ 
            --log-level INFO \
	    VR03.Skera.s_reads.segmented.fltnc.corrected.sorted.dedup.ensemble.mapped.collapsed.gff \
            -o VR03.Skera.s_reads.segmented.fltnc.corrected.sorted.dedup.ensemble.mapped.collapsed.sorted.gff
# pigeon version v1.0.0
# Output: VR03.Skera.s_reads.segmented.fltnc.corrected.sorted.dedup.ensemble.mapped.collapsed.sorted.gff
# Summary output:
# Approximate run time:

#Sort and index the reference files:	
pigeon sort Zea_mays.Zm-B73-REFERENCE-NAM-5.0.56.gtf -o Zea_mays.Zm-B73-REFERENCE-NAM-5.0.56.gtf
pigeon index Zea_mays.Zm-B73-REFERENCE-NAM-5.0.56.gtf
     
#Classify isoforms:
pigeon classify -j 56 \ 
                --log-level INFO \
		VR03.Skera.s_reads.segmented.fltnc.corrected.sorted.dedup.ensemble.mapped.collapsed.sorted.gff \
                Zea_mays.Zm-B73-REFERENCE-NAM-5.0.56.gtf \
		Zea_mays.Zm-B73-REFERENCE-NAM-5.0.dna.toplevel.fa
# Output: VR03_classification.txt; VR03_junctions.txt
# Summary output:
# Approximate run time: 2s (56 threads)

# Filter isoforms:	   
pigeon filter -j 56 \
              --log-level INFO \
	      VR03_classification.txt \
	      -i VR03.Skera.s_reads.segmented.fltnc.corrected.sorted.dedup.ensemble.mapped.collapsed.sorted.gff
# output path: ./ensemble_gft/ 
# *.filtered_lite_classification.txt
# *.filtered_lite_junctions.txt
# *.filtered_lite_reasons.txt
# *.sorted.filtered_lite.gff (only if --isoforms is used)
# Summary output: VR03_classification.filtered.report.json
# Approximate run time: 8s (56 threads)

# Make Seurat compatible input:
pigeon make-seurat -j 56 \
		   --log-level INFO \
                   --dedup VR03.Skera.s_reads.segmented.fltnc.corrected.sorted.dedup.fasta \
		   --group VR03.Skera.s_reads.segmented.fltnc.corrected.sorted.dedup.ensemble.mapped.collapsed.group.txt \
                   -d ensemble_gtf \
		   VR03_classification.filtered_lite_classification.txt
# output path: ./ensembel_gft/ensemble_seurat/
# ./annotated.info.csv
# ./info.csv
# ./genes_seurat/barcodes.tsv
# ./genes_seurat/genes.tsv
# ./genes_seurat/matrix.mtx
# ./isoforms_seurat/barcodes.tsv
# ./isoforms_seurat/genes.tsv
# ./isoforms_seurat/matrix.mtx
# Summary output:
# Approximate run time: 2m 10s (56 threads)
