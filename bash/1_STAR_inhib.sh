#!/bin/bash
#####FASTQ

cd /project/FASTQ

##downsample of fastq files to 100k reads
for i in *R1*.fastq; do name=$(basename ${i} _R1_001.fastq); 
head -q -n 400000 ${name}_R1_001.fastq > /project/trunc/${name}_R1_001_trunc.fastq 
done

#count lines and move files with less than 400k to another different folder (did not pass QC)
wc -l *_R1_001_trunc.fastq



#####STAR alignment

cd /project/trunc/

for i in `ls *R1*.fastq`; do name=$(basename ${i} _R1_001_trunc.fastq); echo STAR --genomeDir /genomes/star/GRCm38.p6.mgp --runMode alignReads --outFilterScoreMinOverLread 0.5 --readFilesIn ${name}_R1_001_trunc.fastq --clip3pAdapterSeq CTGTCTCTTATACACATCT --runThreadN 16 --limitBAMsortRAM 50000000000 --genomeLoad NoSharedMemory  --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM --outReadsUnmapped Fastx --outSAMunmapped Within --outSAMattributes Standard --outFileNamePrefix ./${name}_trunc >> run_STAR.sh;done; 

nohup sh run_STAR.sh >> log_STAR.sh &

tail -F log_STAR.sh


##these were then moved to the /project/STAR


#####index bamfiles
cd /project/STAR
ls *Aligned.sortedByCoord.out.bam | xargs -n1 -P5 samtools index





