#!/bin/bash
#####from bamfile to bw to plotprofle


##making bw bigwig files
cd /project/STAR

for i in *S*Aligned.sortedByCoord.out.bam; do name=$(basename ${i} Aligned.sortedByCoord.out.bam);
bamCoverage -p 15 --normalizeUsing RPGC --effectiveGenomeSize 2652783500 --skipNAs -b ${name}Aligned.sortedByCoord.out.bam -o /project/quant_sort_coord/${name}.bw;
done 



#####deeptools - computematrix
cd /project/quant_sort_coord_bulk/quant_sort_coord;
computeMatrix scale-regions --metagene -R /genomes/mouse/GRCm38.mgp/Mus_musculus.GRCm38.97.chr.gtf -S *.bw -o deeptools/computematrix.mat.gz -p8 --outFileNameMatrix computematrix.tab --verbose



#####deeptools - plotprofile
cd /project/quant_sort_coord_bulk/quant_sort_coord/deeptools;
plotProfile --matrixFile transcriptome_coverage.mat.gz  -out plotprofile.pdf --outFileNameData plotprofile.tab