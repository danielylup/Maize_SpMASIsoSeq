# split the fa sequence of Novel gene model into 100 sequences parts
seqkit split2 -j 56 -s 100 -O ./Maize_V5ext/blastp_analysis/ZmNovelTranscript_sequence2.fa.split/Novel_gene_list/ ./Maize_V5ext/ZmNovelTranscript_sequence.fa

################################################################################################################################################################################
##### remove novel gene models that overlapped with other gene models  #####
#!/usr/bin/env bash

for NUM in {001..207}; do 
    awk '(NR==FNR) { toRemove[$1]; next }
        /^>/ { p=1; for(h in toRemove) if ($0 ~ h) p=0 }
        p' ./Maize_V5ext/gene_list/overlap_NovelGene_list.txt 
        ./Maize_V5ext/blastp_analysis/ZmNovelTranscript_sequence2.fa.split/ZmNovelTranscript_sequence.part_${NUM}.fa > 
        ./Maize_V5ext/blastp_analysis/ZmNovelTranscript_sequence2.fa.split/Novel_gene_list/ZmNovelTranscript2_sequence.part_${NUM}.fa
done

################################################################################################################################################################################
##### Run blastp for all splitted fasta #####
#!/usr/bin/env bash

for NUM in {001..207} do
    ~/tools/ncbi-blast-2.15.0+/bin/blastp -db ~/Project/Reference/BLAST/plant_DB/tmp/nr_plant.fasta 
                                          -query ./Maize_V5ext/blastp_analysis/ZmNovelTranscript_sequence2.fa.split/Novel_gene_list/ZmNovelTranscript2_sequence.part_${NUM}.fa 
                                          -out ./Maize_V5ext/blastp_analysis/blastp_output/plantDB/ZmNovelTranscript_plantDB_blastp_out_${NUM}.txt 
                                          -num_threads 96 -outfmt 6 -max_target_seqs 5
    echo "Finished BLASTP for fasta file ${NUM}."
done

################################################################################################################################################################################
# Calculate number of novel gene model with blastp hit
cat ./Maize_V5ext/blastp_analysis/blastp_output/plantDB/ZmNovelTranscript_plantDB_blastp_out_* | 
  awk '{if($11 <= 1e-5) print}' | 
  cut -f1,2 -d "." | 
  sort | uniq | wc -l



