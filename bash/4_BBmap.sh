#!/bin/bash
####BBMap
#####mhist: Histogram of match, sub, del, and insertion rates

cd /project/trunc/

for i in `ls *S*.fastq`; do name=$(basename ${i} _R1_001_trunc.fastq); echo /programs/bbmap/bbmap.sh -Xmx50g -Xms50g t=24 in=/project/trunc/${name}_R1_001_trunc.fastq  ref=/genomes/Mus_musculus.GRCm38.dna.toplevel.fa mhist=${name}_mhist >> run_bbmap.sh;done;


mv /project/trunc/run_bbmap.sh /project/BBMap

cd /project/BBMap
nohup sh run_bbmap.sh &>> log_bbmap &



#####edit output files for analysis
cd /project/BBMap

mv *_mhist /project/BBMap/mhist

cd /project/BBMap/mhist


for i in *S*_mhist; do name=$(basename ${i} _mhist);
awk '{print $2}' ${name}_mhist | sed "1s/^/${name} \n/" | sed '2d' > ${name}_match;
done

for i in *S*_mhist; do name=$(basename ${i} _mhist);
awk '{print $3}' ${name}_mhist | sed "1s/^/${name} \n/" | sed '2d' > ${name}_sub; 
done

for i in *S*_mhist; do name=$(basename ${i} _mhist);
awk '{print $4}' ${name}_mhist | sed "1s/^/${name} \n/" | sed '2d' > ${name}_del
done

for i in *S*_mhist; do name=$(basename ${i} _mhist);
awk '{print $5}' ${name}_mhist | sed "1s/^/${name} \n/" | sed '2d' > ${name}_ins
done

paste *_match | column -s $'\t' -t > trunc_match_merge
paste *_sub | column -s $'\t' -t > trunc_sub_merge
paste *_del | column -s $'\t' -t > trunc_del_merge
paste *_ins | column -s $'\t' -t > trunc_ins_merge

sed 's/ \+ /\t/g' trunc_match_merge > trunc_match_merge_final
sed 's/ \+ /\t/g' trunc_sub_merge > trunc_sub_merge_final
sed 's/ \+ /\t/g' trunc_del_merge > trunc_del_merge_final
sed 's/ \+ /\t/g' trunc_ins_merge > trunc_ins_merge_final
