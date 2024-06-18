#!/bin/bash
#####Qualimap - RNAseq

cd /project/STAR

for i in `ls *S*Aligned.sortedByCoord.out.bam`; do name=$(basename ${i} Aligned.sortedByCoord.out.bam); echo /programs/qualimap_v2.2.1/qualimap rnaseq -bam /project/STAR/${name}Aligned.sortedByCoord.out.bam -gtf /genomes/Mus_musculus.GRCm38.102.chr.collapse.gtf -oc ${name}_RNA >> run_qualimap_rna.sh;done 

mv /project/STAR/run_qualimap_rna.sh /project/STAR/qualimap_results/

cd /project/STAR/qualimap_results/

nohup sh run_qualimap_rna.sh >> log_qualimap_rna &



#####move "_rnaseq_qc" folders to /qualimap_results/RNAseq/

cd /project/STAR/qualimap_results/

mv -t /project/STAR/qualimap_results/RNAseq /project/STAR/ *Aligned.sortedByCoord.out_rnaseq_qc

cd /project/STAR/qualimap_results/RNAseq



#####edit output file for analysis
find /project/STAR/qualimap_results/RNAseq -type f -name "rnaseq_qc_results.txt" | xargs sed -i '1,13d;23,26d;31,34d;38,41d;43d'

find /project/STAR/qualimap_results/RNAseq -type f -name "rnaseq_qc_results.txt" | xargs sed -i 's/ //g' 

find /project/STAR/qualimap_results/RNAseq -type f -name "rnaseq_qc_results.txt" | xargs sed -i 's/,//g' 

find /project/STAR/qualimap_results/RNAseq -type f -name "rnaseq_qc_results.txt" | xargs sed -i "s/'/_/g" 

find /project/STAR/qualimap_results/RNAseq -type f -name "rnaseq_qc_results.txt" | xargs sed -i "s/=/,/g" 

find /project/STAR/qualimap_results/RNAseq -type f -name "rnaseq_qc_results.txt" | xargs sed -i "s/:/,/g" 

find /project/STAR/qualimap_results/RNAseq -type f -name "rnaseq_qc_results.txt" | xargs sed -i "s/:/,/g" 



##summary 

cd /project/STAR/qualimap_results/RNAseq
 

find . -type f -name "rnaseq_qc_results.txt" -printf "/%P\n" | while read FILE; do DIR=$(dirname "$FILE" );  mv ."$FILE" ."$DIR""$DIR".txt; done 

 
cp -r /project/STAR/qualimap_results/RNAseq/*/*S*.txt /project/STAR/qualimap_results/RNAseq 

cd /project/STAR/qualimap_results/RNAseq/

for f in *Coord.out_rnaseq_qc.txt; do printf '%s\n' 1i "$f" . w q | ed "$f"; done


##remove _truncAligned.sortedByCoord.out_rnaseq_qc.txt inside file (extra bak file)

sed '1 s/_truncAligned.sortedByCoord.out_rnaseq_qc.txt//' *Coord.out_rnaseq_qc.txt -i.bak

##add 'parameters,' to beginning of first line 

sed -i '1s/\(.*\)/parameters,\1/' *Coord.out_rnaseq_qc.txt

sed -i "s/$/,/g" *Coord.out_rnaseq_qc.txt 

sed -i 's/ //g' *Coord.out_rnaseq_qc.txt 

sed -i 's/.$//' *Coord.out_rnaseq_qc.txt 

paste -d"," *Coord.out_rnaseq_qc.txt >> merge_qualimap_RNAseq_trunc1.csv

awk -F',' '{for(i=2;i<=NF;i+=2){printf "%s ",$i;} print ""}' merge_qualimap_RNAseq_trunc1.csv > merge_qualimap_RNAseq_trunc_final1.csv

awk -F',' '{print $1}' merge_qualimap_RNAseq_trunc1.csv > merge_qualimap_RNAseq_trunc_labels1.csv
