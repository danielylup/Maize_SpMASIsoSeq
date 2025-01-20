#HiFi reads processing and isoform structure classification using IsoSeq3 and Pigeon pipeline

skera split -j 56 \
            --log-level INFO \
	    m64531e_230118_204124.hifi_reads.bam \
            mas16_primers.fasta \
	    m64531e.Skera.s_reads.segmented.bam
    #HiFi read segmentation: 
      #skera version 0.1.0  
      #thread number setting available (-j [default 0])
    #output: m64531e.Skera.s_reads.segmented.bam
    #Summary output: m64531e.Skera.s_reads.segmented_2.summary.csv
    #Approximate run time: 4m 23s (56 threads)

lima -j 56 \
     --log-level INFO \
     m64531e.Skera.s_reads.segmented.bam \
     10x_3kit_primers.fasta \
     m64531e.Skera.s_reads.segmented.primer_rm.bam \
     --isoseq
    #Primer removal:
      #lima version 2.6.99 
      #thread number setting available (-j [default 0])
    #output: m64531e.Skera.s_reads.segmented.primer_rm.bam (and others: m64xxxxxx.fl.xxx)		
    #Summary output: m64531e.Skera.s_reads.segmented.primer_rm.lima.summary
    #Approximate run time: 8m 22 (56 threads)

isoseq3 tag -j 56 \
            --log-level INFO \
	    m64531e.Skera.s_reads.segmented.primer_rm.bam \
            m64531e.Skera.s_reads.segmented.flt.bam \
	    --design T-12U-16B
    #Tag clipped:
      #isoseq3 version 3.8.1 
      #thread number setting available (-j [default 0])
    #output: m64531e.Skera.s_reads.segmented.flt.bam (others: m64xxxxxxx.flt.xxx)
    #Summary output: None
    #Approximate run time: 6m 38s (56 threads)

isoseq3 refine -j 56 \
               --log-level INFO \
	       --require-polya \
	       m64531e.Skera.s_reads.segmented.flt.bam \
	       10x_3kit_primers.fasta \
	       m64531e.Skera.s_reads.segmented.fltnc.bam \
    #Refine:    
      #thread number setting available (-j [default 0])
    #output: m64531e.Skera.s_reads.segmented.fltnc.bam
    #Summary output: m64531e.Skera.s_reads.segmented.fltnc.report.csv
    #Approximate run time: 3m 34s (56 threads)
	   
isoseq3 correct -j 56 \
                --log-level INFO \
		--barcodes visium_v1.RC.txt \
                m64531e.Skera.s_reads.segmented.fltnc.bam \
		m64531e.Skera.s_reads.segmented.fltnc.corrected.bam
    #Cell barcode correction: 
      #thread number setting available (-j [default 0])
    #output: m64531e.Skera.s_reads.segmented.fltnc.corrected.bam
    #Summary output: m64531e.Skera.s_reads.segmented.fltnc.corrected.report.json
    #Approximate run time: 335s (56 threads)
	   
samtools sort -@ 56 \
              -t CB m64531e.Skera.s_reads.segmented.fltnc.corrected.bam \
	      -o m64531e.Skera.s_reads.segmented.fltnc.corrected.sorted.bam
isoseq3 groupdedup -j 56 \
                   --log-level INFO \ 
		   m64531e.Skera.s_reads.segmented.fltnc.corrected.sorted.bam \
                   m64531e.Skera.s_reads.segmented.fltnc.corrected.sorted.dedup.bam
    #Deduplication: 
      #thread number setting available (-@ & -j [default 0])
    #output: m64531e.Skera.s_reads.segmented.fltnc.corrected.sorted.dedup.bam
    #Summary output: m64531e.Skera.s_reads.segmented.fltnc.corrected.report.json
    #Approximate run time: 22m 27s (56 threads; not including samtools)
	   
pbmm2 align -j 56 \
            --log-level INFO \
	    --preset ISOSEQ \
            --sort m64531e.Skera.s_reads.segmented.fltnc.corrected.sorted.dedup.bam \
	    Zea_mays.Zm-B73-REFERENCE-NAM-5.0.dna.toplevel.fa \
            m64531e.Skera.s_reads.segmented.fltnc.corrected.sorted.dedup.ensemble.mapped.bam 
    #Map reads to a reference genome:
      #pbmm2 version 1.9.0
    #output: m64531e.Skera.s_reads.segmented.fltnc.corrected.sorted.dedup.ensemble.mapped.bam
    #Summary output:
    #Approximate run time:
	   
isoseq3 collapse -j 56 \
		 --log-level INFO \
                 m64531e.Skera.s_reads.segmented.fltnc.corrected.sorted.dedup.ensemble.mapped.bam \
		 m64531e.Skera.s_reads.segmented.fltnc.corrected.sorted.dedup.ensemble.mapped.collapsed.gff
    #Collapse into unique isoforms: 
    #output: m64531e.Skera.s_reads.segmented.fltnc.corrected.sorted.dedup.ensemble.mapped.collapsed.gff
    #Summary output: 
    #Approximate run time: 3m 27s
	   
pigeon sort -j 56 \ 
            --log-level INFO \
	    m64531e.Skera.s_reads.segmented.fltnc.corrected.sorted.dedup.ensemble.mapped.collapsed.gff \
            -o m64531e.Skera.s_reads.segmented.fltnc.corrected.sorted.dedup.ensemble.mapped.collapsed.sorted.gff
    #Sort input transcript GFF: 
      #pigeon version v1.0.0
    #output: m64531e.Skera.s_reads.segmented.fltnc.corrected.sorted.dedup.ensemble.mapped.collapsed.sorted.gff
    #Summary output:
    #Approximate run time:
	
pigeon sort Zea_mays.Zm-B73-REFERENCE-NAM-5.0.56.gtf -o Zea_mays.Zm-B73-REFERENCE-NAM-5.0.56.gtf
pigeon index Zea_mays.Zm-B73-REFERENCE-NAM-5.0.56.gtf
    #Sort and index the reference files: 

pigeon classify -j 56 \ 
                --log-level INFO \
		m64531e.Skera.s_reads.segmented.fltnc.corrected.sorted.dedup.ensemble.mapped.collapsed.sorted.gff \
                Zea_mays.Zm-B73-REFERENCE-NAM-5.0.56.gtf \
		Zea_mays.Zm-B73-REFERENCE-NAM-5.0.dna.toplevel.fa
    #Classify isoforms: 
    #output: m64531e_classification.txt; m64531e_junctions.txt
    #Summary output:
    #Approximate run time: 2s (56 threads)
	   
pigeon filter -j 56 \
              --log-level INFO \
	      m64531e_classification.txt \
	      -i m64531e.Skera.s_reads.segmented.fltnc.corrected.sorted.dedup.ensemble.mapped.collapsed.sorted.gff
    #Filter isoforms: 
    #output path: ./ensemble_gft/ 
    	#*.filtered_lite_classification.txt
        #*.filtered_lite_junctions.txt
        #*.filtered_lite_reasons.txt
        #*.sorted.filtered_lite.gff (only if --isoforms is used)
    #Summary output: m64531e_classification.filtered.report.json
    #Approximate run time: 8s (56 threads)

pigeon make-seurat -j 56 \
		   --log-level INFO \
                   --dedup m64531e.Skera.s_reads.segmented.fltnc.corrected.sorted.dedup.fasta \
		   --group m64531e.Skera.s_reads.segmented.fltnc.corrected.sorted.dedup.ensemble.mapped.collapsed.group.txt \
                   -d ensemble_gtf \
		   m64531e_classification.filtered_lite_classification.txt
    #Make Seurat compatible input:
    #output path: ./ensembel_gft/ensemble_seurat/
        #./annotated.info.csv
        #./info.csv
        #./genes_seurat/barcodes.tsv
        #./genes_seurat/genes.tsv
        #./genes_seurat/matrix.mtx
        #./isoforms_seurat/barcodes.tsv
        #./isoforms_seurat/genes.tsv
        #./isoforms_seurat/matrix.mtx
    #Summary output:
    #Approximate run time: 2m 10s (56 threads)
