# Using the python script provided in CPAT tool, the nucleotide sequences of novel gene model were translated to amino acid sequences 
cpat.py -x Zm-B73_Hexamer.tsv --antisense -d Zm-B73.logit.RData --top-orf=5 -g NovelGene_filtered_lite.fa -o NovelGene
# the default setting of "best" ORF selection is based on the probabilty evaluation of several criteria, including Fickett Score, Hexamer Score, and length.
# the model used for the evalution (Zm-B73_Hexamer.tsv) is trained from the coding and non-coding sequencing files obtained from NCBI (https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/902/167/145/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0/) 

###################################################################################################################################################################################################################################################
###### Run PFAMscan on the novel gene AA sequence files ######
#!/usr/bin/env bash
set -e

for NUM in {001..207}; do
        pfam_scan.pl -fasta ./Maize_V5ext/blastp_analysis/ZmNovelTranscript_sequence_AA_split/NovelGene_best_ORF_AA.part_${NUM}.fa 
                     -outfile ./Maize_V5ext/PFAMscan_analysis/PFAMscan_output/ZmNovelTranscript_pfamscan_out_${NUM}.txt 
                     --dir ./ncRNA_analysis/CPAT_analysis/PfamA_HMM/
        echo "Finished pfam scanning for sequence file ${NUM}!"
done

###################################################################################################################################################################################################################################################
### Calculate how many novel gene model with PFAM support
grep -v "#" ./Maize_V5ext/PFAMscan_analysis/PFAMscan_output/ZmNovelTranscript_pfamscan_out_*.txt | cut -f2 -d ":" | sed '/^[[:space:]]*$/d' | cut -f1,2 -d"_" | sort | uniq | wc -l
